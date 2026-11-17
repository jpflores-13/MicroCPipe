## Shared helpers for MicroCPipe.snakefile and MicroCMerge.snakefile.
## Rules carry no cpu/memory/time of their own beyond `threads`; resources live in
## profiles/<site>/config.yaml (set-threads / set-resources), as in dietJuicer.

import os
import sys

OUTPUTS = list(config.get("outputs", ["cool", "mcool", "hic"]))
_bad = set(OUTPUTS) - {"cool", "mcool", "hic"}
if _bad:
    raise ValueError(f"config 'outputs' has unknown format(s) {sorted(_bad)}; use cool, mcool, hic")

BASE_RES = int(config.get("base_resolution", 1000))
RESOLUTIONS = sorted({int(r) for r in config.get("resolutions", [BASE_RES])} | {BASE_RES})
for _r in RESOLUTIONS:
    if _r % BASE_RES:
        raise ValueError(f"resolution {_r} is not a multiple of base_resolution {BASE_RES}")
HIC_RES = sorted({int(r) for r in config.get("hic_resolutions", RESOLUTIONS)})


def outdir():
    """outputDir, with a relative path placed under $MICROC_OUTPUT_ROOT when set."""
    d = os.path.expanduser(str(config.get("outputDir", "output")))
    if not os.path.isabs(d):
        d = os.path.join(os.environ.get("MICROC_OUTPUT_ROOT", os.getcwd()), d)
    return d.rstrip("/")


OUT = outdir()
ENVS = os.path.join(workflow.basedir, "..", "envs")


def map_targets(groups):
    """The contact maps requested in config['outputs'] (what QC and benchmarking wait for)."""
    t = []
    for g in groups:
        if "cool" in OUTPUTS:
            t.append(f"{OUT}/{g}/{g}.cool")
        if "mcool" in OUTPUTS:
            t.append(f"{OUT}/{g}/{g}.mcool")
        if "hic" in OUTPUTS:
            t.append(f"{OUT}/{g}/{g}.hic")
    return t


def analysis_targets(groups):
    """Optional downstream analyses requested in config['analysis'] (never block QC)."""
    return [f"{OUT}/{g}/analysis/{g}_{a}.done" for g in groups for a in config.get("analysis", [])]


def format_targets(groups):
    """Everything `rule all` asks for: the maps plus any optional analyses."""
    return map_targets(groups) + analysis_targets(groups)


##### Email notifications (see `notify:` in the config) #####
NOTIFY = config.get("notify", {}) or {}
EMAIL = os.environ.get("MICROC_EMAIL") or NOTIFY.get("email") or ""


def notify(event, body, run_name):
    """Send one message for a workflow event, if an address and the event are configured."""
    if not EMAIL or event not in NOTIFY.get("events", ["start", "success", "error"]):
        return
    import shlex, socket, subprocess
    subject = f"{run_name}: {event}"
    body = (f"Run:     {run_name}\nOutput:  {OUT}\nHost:    {socket.gethostname()}\n\n{body}\n")
    cmd = NOTIFY.get("command", "mail -s {subject} {email}").format(
        subject=shlex.quote(subject), email=shlex.quote(EMAIL))
    try:
        subprocess.run(cmd, shell=True, input=body, text=True, timeout=60, check=True)
    except Exception as e:  # a failed email must never fail the workflow
        print(f"Could not send '{event}' email to {EMAIL}: {e}", file=sys.stderr)


def log_tail(path, n=40):
    """Last n lines of the Snakemake log, for the failure email."""
    try:
        with open(path) as fh:
            return "".join(fh.readlines()[-n:])
    except (OSError, TypeError):
        return "(Snakemake log unavailable)"
