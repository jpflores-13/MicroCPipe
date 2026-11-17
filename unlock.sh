#! /bin/bash -login

## Release the Snakemake lock after a failed or interrupted run, then
## resubmit with MicroCPipe.sh (or MicroCMerge.sh).
##
## Usage:
##   bash unlock.sh                                  # core workflow, default config
##   bash unlock.sh --merge                          # merge workflow
##   bash unlock.sh --configfile other.yaml --profile profiles/<site>
## The profile is chosen as in MicroCPipe.sh: --profile, $MICROC_PROFILE,
## then hpcProfile in the config.

## Exit if any command fails
set -e

CONFIG="config/MicroCConfig.yaml"
PROFILE=""
SNAKEFILE="workflows/MicroCPipe.snakefile"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --configfile) CONFIG="$2";  shift 2 ;;
        --profile)    PROFILE="$2"; shift 2 ;;
        --merge)      SNAKEFILE="workflows/MicroCMerge.snakefile"; shift ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

cfg() { sed -nE "s/^$1:[[:space:]]*['\"]?([^'\"#[:space:]]+).*/\1/p" "${CONFIG}" | head -1; }
PROFILE="${PROFILE:-${MICROC_PROFILE:-$(cfg hpcProfile)}}"

## Load this site's software
[[ -f "${PROFILE%/}/env.sh" ]] || { echo "No ${PROFILE%/}/env.sh: set --profile, $MICROC_PROFILE or hpcProfile in ${CONFIG}" >&2; exit 1; }
source "${PROFILE%/}/env.sh"

## Unlock snakemake workflow (venv is not activated; see MicroCPipe.sh)
env/bin/python3 -m snakemake --unlock -s "${SNAKEFILE}" --cores 1 --configfile "${CONFIG}"

## Success message
echo -e "\033[0;32mDirectory unlocked, ready to rerun.\033[0m"
