#!/usr/bin/env bash
# to_pdf.sh — render pipeline HTML reports to PDF. Uses whatever PDF engine is available
# (headless Chrome/Chromium, wkhtmltopdf, or weasyprint). Pass HTML files or a results dir.
#
#   bin/to_pdf.sh results/*/*.report.html results/*.cohort.html      # explicit files
#   bin/to_pdf.sh results                                            # whole results tree
set -euo pipefail

# collect inputs
inputs=()
for a in "$@"; do
  if [ -d "$a" ]; then
    while IFS= read -r f; do inputs+=("$f"); done < <(find "$a" -name '*.report.html' -o -name '*.cohort.html')
  else inputs+=("$a"); fi
done
[ ${#inputs[@]} -eq 0 ] && { echo "usage: to_pdf.sh <html files | results dir>"; exit 1; }

# pick an engine
CHROME=""
for c in "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
         "/Applications/Chromium.app/Contents/MacOS/Chromium" \
         "$(command -v google-chrome 2>/dev/null||true)" "$(command -v chromium 2>/dev/null||true)" \
         "$(command -v chromium-browser 2>/dev/null||true)"; do
  [ -n "$c" ] && [ -x "$c" ] && CHROME="$c" && break
done
WK="$(command -v wkhtmltopdf 2>/dev/null || true)"
WE="$(command -v weasyprint 2>/dev/null || true)"

render() {  # $1 html -> ${1%.html}.pdf
  local html="$1" pdf="${1%.html}.pdf"
  if [ -n "$CHROME" ]; then
    "$CHROME" --headless --disable-gpu --no-pdf-header-footer --print-to-pdf="$pdf" "file://$(cd "$(dirname "$html")"&&pwd)/$(basename "$html")" 2>/dev/null
  elif [ -n "$WK" ]; then "$WK" -q "$html" "$pdf"
  elif [ -n "$WE" ]; then "$WE" "$html" "$pdf"
  else echo "no PDF engine (install Chrome, wkhtmltopdf, or weasyprint)"; exit 2; fi
  echo "  -> $pdf"
}

echo "[to_pdf] engine: ${CHROME:-${WK:-${WE:-none}}}"
for f in "${inputs[@]}"; do [ -f "$f" ] && render "$f"; done
