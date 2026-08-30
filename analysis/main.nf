#!/usr/bin/env nextflow
/*
 * MalaCall — end-to-end pan-Plasmodium analysis & reporting pipeline (Nextflow DSL2)
 *
 *   MAP → SPECIES_ID → { RESISTANCE, CNV_DELETION, DIVERSITY, PK_CLUSTER } → REPORT → COHORT
 *
 * The heavy lifting is in the Python scripts under bin/ (auto-added to PATH by Nextflow); this file is
 * the orchestration + provenance layer. Every stage is validated against simulated ground truth
 * (see test/validate.py). Run:  nextflow run main.nf --samplesheet samples.csv --outdir results
 */
nextflow.enable.dsl = 2

params.samplesheet = null
params.outdir      = 'results'
params.ref         = "${projectDir}/../refs/pan_reference.fa.gz"
params.manifest    = "${projectDir}/manifest/amplicons.tsv"
params.codon_map   = "${projectDir}/../targets/codon_to_genomic_pan.tsv"
params.pk_panel    = "${projectDir}/refs/pk_diag_panel.tsv"
params.platform    = 'ont'          // 'ont' | 'illumina'
params.cohort_name = 'cohort'
// --- custom-barcode demultiplexing (ONT/PCR combinatorial 24-nt barcodes; MinKNOW can't do these) ---
params.demux       = false          // set true to demultiplex a pooled run before analysis
params.pooled      = null           // pooled FASTQ(s) (glob) when --demux
params.barcodes    = "${projectDir}/../refs/barcodes/barcodes_96.tsv"
params.barcode_map = null           // CSV: sample,fwd_bc,rev_bc
params.demux_mode  = 'dual'         // 'dual' | 'single'
params.bc_mismatch = 2
params.python      = 'python3'      // interpreter with pysam/numpy/pandas (e.g. the captureforge env)
params.min_frac    = 0.01
params.min_vaf     = 0.05

def BIN = "${projectDir}/bin"
def PY  = { "PYTHONPATH=${BIN} ${params.python}" }

process DEMUX {
    publishDir "${params.outdir}/demux", mode: 'copy'
    input:  path(pooled); path(barcodes); path(bcmap)
    output: path("demuxed/*.fastq.gz"), emit: reads
            path("demuxed/demux_stats.tsv")
    script:
    """
    ${PY()} ${BIN}/demux.py --fastq ${pooled} --barcodes ${barcodes} --sample-map ${bcmap} \
        --outdir demuxed --mode ${params.demux_mode} --bc-mismatch ${params.bc_mismatch} --write-samplesheet
    rm -f demuxed/unclassified.fastq.gz            # the unclassified bin is not analysed
    """
}

process MAP {
    tag "$sample"
    publishDir "${params.outdir}/${sample}", mode: 'copy', pattern: '*.bam*'
    input:  tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")
    script:
    def preset = params.platform == 'ont' ? 'map-ont' : 'sr'
    def reads  = r2.name != 'NO_FILE' ? "${r1} ${r2}" : "${r1}"
    """
    minimap2 -ax ${preset} -t ${task.cpus} --secondary=no ${params.ref} ${reads} \
      | samtools sort -@ ${task.cpus} -o ${sample}.bam
    samtools index ${sample}.bam
    """
}

process SPECIES_ID {
    tag "$sample"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    input:  tuple val(sample), path(bam), path(bai)
    output: tuple val(sample), path(bam), path(bai), path("${sample}.species.tsv"), path("${sample}.species.call")
    script:
    """
    ${PY()} ${BIN}/species_id.py ${bam} ${params.manifest} ${sample} ${sample}.species.tsv --min-frac ${params.min_frac}
    """
}

process RESISTANCE {
    tag "$sample"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    input:  tuple val(sample), path(bam), path(bai), path(species_tsv), path(call)
    output: tuple val(sample), path("${sample}.resistance.tsv"), path("${sample}.drugs.tsv")
    script:
    """
    SP=\$(tr '\\t' ',' < ${call}); [ -z "\$SP" ] && SP="Pf,Pv,Poc,Pow,Pm,Pk"
    ${PY()} ${BIN}/resistance_genotype.py ${bam} ${params.manifest} ${params.codon_map} ${sample} . \
        --species "\$SP" --min-vaf ${params.min_vaf}
    """
}

process CNV_DELETION {
    tag "$sample"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    input:  tuple val(sample), path(bam), path(bai), path(species_tsv), path(call)
    output: tuple val(sample), path("${sample}.cnv.tsv")
    script:
    """
    SP=\$(tr '\\t' ',' < ${call}); [ -z "\$SP" ] && SP="Pf,Pv,Poc,Pow,Pm,Pk"
    ${PY()} ${BIN}/cnv_deletion.py ${bam} ${params.manifest} ${sample} ${sample}.cnv.tsv --species "\$SP"
    """
}

process DIVERSITY {
    tag "$sample"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    input:  tuple val(sample), path(bam), path(bai), path(species_tsv), path(call)
    output: tuple val(sample), path("${sample}.diversity.tsv"), path("${sample}.microhap.tsv"), path("${sample}.alleles.tsv"), path("${sample}.consensus.tsv")
    script:
    """
    SP=\$(tr '\\t' ',' < ${call}); [ -z "\$SP" ] && SP="Pf,Pv,Poc,Pow,Pm,Pk"
    ${PY()} ${BIN}/diversity.py ${bam} ${params.manifest} ${sample} . --species "\$SP"
    """
}

process PK_CLUSTER {
    tag "$sample"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    input:  tuple val(sample), path(bam), path(bai), path(species_tsv), path(call)
    output: tuple val(sample), path("${sample}.cluster.tsv"), optional: true
    script:
    // runs for every sample; the script writes nothing unless P. knowlesi reads are present
    """
    if grep -qx Pk <(tr '\\t' '\\n' < ${call}); then
        ${PY()} ${BIN}/pk_cluster.py ${bam} ${params.manifest} ${params.pk_panel} ${sample} ${sample}.cluster.tsv
    fi
    """
}

process REPORT {
    tag "$sample"
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    input:  tuple val(sample), path('*')
    output: path("${sample}.report.html")
    script:
    """
    ${PY()} ${BIN}/report.py . ${sample}
    """
}

process COHORT {
    publishDir "${params.outdir}", mode: 'copy'
    input:  path(files, stageAs: 'in/*')
    output: path("${params.cohort_name}.cohort.*")
    script:
    """
    mkdir -p cohort
    for f in in/*; do s=\$(basename \$f | cut -d. -f1); mkdir -p cohort/\$s; cp \$f cohort/\$s/; done
    ${PY()} ${BIN}/cohort_report.py cohort ${params.cohort_name}
    cp cohort/${params.cohort_name}.cohort.* .
    """
}

process RELATEDNESS {
    publishDir "${params.outdir}", mode: 'copy'
    input:  path(alleles, stageAs: 'in/*')
    output: path("${params.cohort_name}.coi.tsv"); path("${params.cohort_name}.relatedness.tsv")
            path("${params.cohort_name}.related_pairs.tsv"); path("export/*")
    script:
    // COI-MLE + pairwise relatedness across the cohort + exports for THE REAL McCOIL / Dcifer
    """
    mkdir -p cohort
    for f in in/*; do s=\$(basename \$f | cut -d. -f1); mkdir -p cohort/\$s; cp \$f cohort/\$s/; done
    ${PY()} ${BIN}/relatedness.py cohort ${params.cohort_name}
    cp cohort/${params.cohort_name}.coi.tsv cohort/${params.cohort_name}.relatedness.tsv cohort/${params.cohort_name}.related_pairs.tsv .
    cp -r cohort/export .
    """
}

workflow {
    if (params.demux) {
        // pooled run -> demultiplex our custom 24-nt barcodes -> per-sample reads
        DEMUX(Channel.fromPath(params.pooled).collect(), file(params.barcodes), file(params.barcode_map))
        reads_ch = DEMUX.out.reads.flatten()
            .map { f -> tuple(f.getBaseName().replaceAll(/\.fastq$/, ''), f, file('NO_FILE')) }
    } else {
        reads_ch = Channel.fromPath(params.samplesheet)
            | splitCsv(header: true)
            | map { row -> tuple(row.sample, file(row.fastq), row.fastq2 ? file(row.fastq2) : file('NO_FILE')) }
    }

    MAP(reads_ch)
    SPECIES_ID(MAP.out)
    RESISTANCE(SPECIES_ID.out)
    CNV_DELETION(SPECIES_ID.out)
    DIVERSITY(SPECIES_ID.out)
    PK_CLUSTER(SPECIES_ID.out)

    // gather every per-sample output for the per-sample report
    per_sample = SPECIES_ID.out.map { s, b, i, tsv, c -> tuple(s, [tsv]) }
        .join(RESISTANCE.out.map    { s, r, d -> tuple(s, [r, d]) })
        .join(CNV_DELETION.out.map  { s, c -> tuple(s, [c]) })
        .join(DIVERSITY.out.map     { s, d, m, al, co -> tuple(s, [d, m]) })
        .join(PK_CLUSTER.out.map    { s, c -> tuple(s, [c]) }, remainder: true)
        .map { row -> tuple(row[0], (row[1] + row[2] + row[3] + row[4] + (row[5] ?: [])).flatten()) }
    REPORT(per_sample)

    // cohort COI-MLE + pairwise relatedness (+ McCOIL / Dcifer exports) over all microhaplotype calls
    RELATEDNESS(DIVERSITY.out.flatMap { s, d, m, al, co -> [al, co] }.collect())

    // cohort surveillance summary over all per-sample tables
    all_tables = SPECIES_ID.out.map   { s, b, i, t, c -> t }
        .mix(RESISTANCE.out.flatMap   { s, r, d -> [r, d] })
        .mix(CNV_DELETION.out.map     { s, c -> c })
        .mix(DIVERSITY.out.map        { s, d, m -> d })
        .collect()
    COHORT(all_tables)
}
