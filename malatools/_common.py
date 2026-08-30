"""Shared helpers for the malaforge / malacall CLIs: locate the repo, run bundled Python
scripts with the current interpreter, and drive Nextflow with a sane environment."""
import os, sys, shutil, subprocess

# malatools/ lives at the repo root, so the repo is one directory up from this file.
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

C_RED = "\033[38;5;131m"; C_TEAL = "\033[38;5;30m"; C_DIM = "\033[2m"; C_OFF = "\033[0m"
def _c(s, col):
    return f"{col}{s}{C_OFF}" if sys.stdout.isatty() else s


def script(relpath):
    """Absolute path to a bundled script, checked to exist."""
    p = os.path.join(REPO, relpath)
    if not os.path.exists(p):
        sys.exit(f"malatools: bundled script not found: {relpath} (expected under {REPO})")
    return p


def run_script(relpath, argv):
    """Run a bundled Python script with the SAME interpreter that runs this CLI (so pysam/
    numpy/pandas are available), forwarding argv verbatim. Returns the child exit code."""
    cmd = [sys.executable, script(relpath), *argv]
    return subprocess.call(cmd)


def find_nextflow(explicit=None):
    for cand in (explicit, os.environ.get("MALACALL_NEXTFLOW"),
                 shutil.which("nextflow"), os.path.expanduser("~/nextflow")):
        if cand and os.path.exists(os.path.expanduser(cand)):
            return os.path.expanduser(cand)
    return None


def java_home(explicit=None):
    # Explicit wins; then a deliberate MALACALL_JAVA_HOME. Prefer the real system JDK from the macOS
    # helper OVER an inherited JAVA_HOME, because conda envs inject their own (often Nextflow-incompatible)
    # JVM into JAVA_HOME. Fall back to JAVA_HOME last (e.g. on Linux where the helper is absent).
    for cand in (explicit, os.environ.get("MALACALL_JAVA_HOME")):
        if cand and os.path.isdir(cand):
            return cand
    try:
        out = subprocess.check_output(["/usr/libexec/java_home"], text=True,
                                      stderr=subprocess.DEVNULL).strip()
        if out and os.path.isdir(out):
            return out
    except Exception:
        pass
    jh = os.environ.get("JAVA_HOME")
    return jh if jh and os.path.isdir(jh) else None


def top_help(prog, tagline, commands, epilog=""):
    """Render a clean, aligned top-level help screen (command → description)."""
    lines = [_c(prog, C_TEAL) + " — " + tagline, ""]
    lines.append(_c("usage:", C_DIM) + f" {prog} <command> [options]")
    lines.append("")
    lines.append(_c("commands:", C_DIM))
    w = max(len(k) for k, _ in commands)
    for k, desc in commands:
        lines.append(f"  {_c(k.ljust(w), C_TEAL)}  {desc}")
    lines.append("")
    lines.append(f"Run '{prog} <command> --help' for options on a specific command.")
    lines.append(f"'{prog} --version' prints the version.")
    if epilog:
        lines += ["", _c(epilog, C_DIM)]
    return "\n".join(lines)


def dispatch(prog, tagline, commands, handlers, epilog=""):
    """Generic sub-command dispatcher.
      commands : ordered list of (name, description) for help
      handlers : {name: relpath_str  OR  callable(argv)->int}
    """
    argv = sys.argv[1:]
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(top_help(prog, tagline, commands, epilog)); return 0
    if argv[0] in ("-V", "--version", "version"):
        from . import __version__
        print(f"{prog} {__version__}"); return 0
    cmd, rest = argv[0], argv[1:]
    h = handlers.get(cmd)
    if h is None:
        print(top_help(prog, tagline, commands, epilog))
        sys.exit(f"\n{prog}: unknown command '{cmd}'")
    if callable(h):
        return h(rest) or 0
    return run_script(h, rest)     # h is a script relpath → forward args to it
