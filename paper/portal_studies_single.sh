#!/usr/bin/env bash
# paper/portal_studies_single.sh
#
# Submit the portal-study benchmark as a SINGLE bsub job (sequential
# walk through every YAML in --configs-dir, or every study folder under
# --studies-root). Useful when:
#   * the array-mode dispatcher (portal_studies_hpc.sh) is overkill;
#   * you want to confine all studies to one host (cache locality);
#   * you want a single bsub log file for the whole sweep.
#
# Memory / wall / cores / queue are all configurable; LSF's `#BSUB`
# directives in portal_studies.bsub are overridden via bsub flags.
#
# Usage:
#   ./paper/portal_studies_single.sh --configs-dir paper/configs
#   ./paper/portal_studies_single.sh --configs-dir <path> --mem 64000 --wall 12:00
#   ./paper/portal_studies_single.sh --configs-dir <path> --only 2327,2326
#   ./paper/portal_studies_single.sh --studies-root /research/.../studies
#
# Flags (all optional except one of --configs-dir / --studies-root):
#   --configs-dir <dir>   Folder of YAML configs.
#   --studies-root <dir>  Folder of per-study subfolders (legacy mode).
#   --only a,b,c          CSV of study IDs to restrict to.
#   --mem <MB>            Memory request (default: 32000).
#   --cores <n>           Cores per job (default: 4).
#   --wall <hh:mm>        Wall-clock limit (default: 24:00).
#   --queue <name>        LSF queue (default: standard).
#   --project <id>        Charge code (default: scminer).
#   --dry-run             Print the bsub command without submitting.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CONFIGS_DIR=""
STUDIES_ROOT=""
ONLY=""
MEM=32000
CORES=4
WALL="24:00"
QUEUE="standard"
PROJECT="scminer"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --configs-dir)  CONFIGS_DIR="$2"; shift 2 ;;
    --studies-root) STUDIES_ROOT="$2"; shift 2 ;;
    --only)         ONLY="$2"; shift 2 ;;
    --mem)          MEM="$2"; shift 2 ;;
    --cores)        CORES="$2"; shift 2 ;;
    --wall|-W)      WALL="$2"; shift 2 ;;
    --queue)        QUEUE="$2"; shift 2 ;;
    --project)      PROJECT="$2"; shift 2 ;;
    --dry-run)      DRY_RUN=1; shift ;;
    -h|--help)
      sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "unknown flag: $1" >&2; exit 2 ;;
  esac
done

if [[ -n "$CONFIGS_DIR" && -n "$STUDIES_ROOT" ]]; then
  echo "ERROR: pass --configs-dir OR --studies-root, not both" >&2
  exit 2
fi
if [[ -z "$CONFIGS_DIR" && -z "$STUDIES_ROOT" ]]; then
  echo "ERROR: pass --configs-dir <dir> or --studies-root <dir>" >&2
  exit 2
fi

mkdir -p paper/metrics paper/logs

STAMP="$(date +%Y%m%d_%H%M%S)"

# Build the env string handed to the .bsub task. -env "all, K=V, K=V"
# starts from the submitter's env, then sets the listed keys.
ENV_PARTS=("all")
if [[ -n "$CONFIGS_DIR" ]]; then
  CONFIGS_ABS="$(cd "$CONFIGS_DIR" && pwd)"
  ENV_PARTS+=("CONFIGS_DIR=${CONFIGS_ABS}")
fi
if [[ -n "$STUDIES_ROOT" ]]; then
  STUDIES_ABS="$(cd "$STUDIES_ROOT" && pwd)"
  ENV_PARTS+=("STUDIES_ROOT=${STUDIES_ABS}")
fi
if [[ -n "$ONLY" ]]; then
  ENV_PARTS+=("ONLY=${ONLY}")
fi
ENV_STR="$(IFS=,; echo "${ENV_PARTS[*]}")"

BSUB_CMD=(
  bsub
    -J "scminer_portal_single_${STAMP}"
    -P "$PROJECT"
    -q "$QUEUE"
    -n "$CORES"
    -R "rusage[mem=${MEM}]"
    -M "${MEM}"
    -W "$WALL"
    -o "paper/logs/portal_single_${STAMP}.out"
    -e "paper/logs/portal_single_${STAMP}.err"
    -env "$ENV_STR"
)

echo "Submitting single bsub job (sequential walk)"
echo "  queue:       $QUEUE"
echo "  cores:       $CORES"
echo "  mem (MB):    $MEM"
echo "  wall:        $WALL"
if [[ -n "$CONFIGS_DIR" ]]; then echo "  configs-dir: $CONFIGS_ABS"; fi
if [[ -n "$STUDIES_ROOT" ]]; then echo "  studies-root: $STUDIES_ABS"; fi
if [[ -n "$ONLY" ]]; then echo "  only:        $ONLY"; fi
echo

if [[ "$DRY_RUN" == "1" ]]; then
  echo "[dry-run]:"
  printf '    %q ' "${BSUB_CMD[@]}"; echo "< paper/portal_studies.bsub"
  exit 0
fi

"${BSUB_CMD[@]}" < paper/portal_studies.bsub

echo
echo "Watch progress:"
echo "    bjobs -J scminer_portal_single_${STAMP}"
echo "    tail -f paper/logs/portal_single_${STAMP}.out"
echo "Final TSV: paper/metrics/portal_studies.tsv"
