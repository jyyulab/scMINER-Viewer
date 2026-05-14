#!/usr/bin/env bash
# paper/portal_studies_hpc.sh
#
# One-command driver for the portal-study benchmark on St. Jude HPC.
#
# Submits an LSF job array (one task per study) so all 21 studies run
# in parallel, then submits a dependent "merge" job that concatenates
# the per-study TSVs into paper/metrics/portal_studies.tsv once the
# array has finished.
#
# Two input layouts are supported, pick exactly one:
#
#   --configs-dir <dir>   folder of N YAML configs (preferred)
#   --studies-root <dir>  folder of N per-study subfolders
#
# Usage:
#   ./paper/portal_studies_hpc.sh --configs-dir paper/configs
#   ./paper/portal_studies_hpc.sh --studies-root /research/.../studies
#   ./paper/portal_studies_hpc.sh --configs-dir <path> --mem 64000 -W 12:00
#
# Flags:
#   --configs-dir <dir>    Folder of YAML configs (one per study).
#   --studies-root <dir>   Folder of per-study subfolders (legacy mode).
#   --queue <name>         LSF queue (default: standard).
#   --mem <MB>             Memory request per task (default: 32000).
#   --cores <n>            Cores per task (default: 4).
#   --wall <hh:mm>         Wall-clock limit per task (default: 6:00).
#   --project <id>         Charge code (default: scminer).
#   --dry-run              Print bsub commands without submitting.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CONFIGS_DIR=""
STUDIES_ROOT=""
QUEUE="standard"
MEM=32000
CORES=4
WALL="6:00"
PROJECT="scminer"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --configs-dir)  CONFIGS_DIR="$2"; shift 2 ;;
    --studies-root) STUDIES_ROOT="$2"; shift 2 ;;
    --queue)        QUEUE="$2"; shift 2 ;;
    --mem)          MEM="$2"; shift 2 ;;
    --cores)        CORES="$2"; shift 2 ;;
    --wall|-W)      WALL="$2"; shift 2 ;;
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

mkdir -p paper/metrics paper/logs paper/hpc

STAMP="$(date +%Y%m%d_%H%M%S)"
MANIFEST_FILE="paper/hpc/manifest_${STAMP}.txt"

if [[ -n "$CONFIGS_DIR" ]]; then
  if [[ ! -d "$CONFIGS_DIR" ]]; then
    echo "ERROR: configs-dir does not exist: $CONFIGS_DIR" >&2
    exit 2
  fi
  # One YAML config path per line, absolute, sorted.
  CONFIGS_ABS="$(cd "$CONFIGS_DIR" && pwd)"
  find "$CONFIGS_ABS" -mindepth 1 -maxdepth 1 -type f \
       \( -name '*.yaml' -o -name '*.yml' \) \
    | sort > "$MANIFEST_FILE"
  N=$(wc -l < "$MANIFEST_FILE" | tr -d ' ')
  if [[ "$N" -eq 0 ]]; then
    echo "ERROR: no *.yaml in $CONFIGS_ABS" >&2
    exit 1
  fi
  echo "Found $N YAML configs in $CONFIGS_ABS"
else
  if [[ ! -d "$STUDIES_ROOT" ]]; then
    echo "ERROR: studies-root does not exist: $STUDIES_ROOT" >&2
    exit 2
  fi
  # One study folder per line (must contain config.yaml or expression.rds).
  STUDIES_ABS="$(cd "$STUDIES_ROOT" && pwd)"
  find "$STUDIES_ABS" -mindepth 1 -maxdepth 1 -type d \
    | sort \
    | while read -r d; do
        if [[ -f "$d/config.yaml" || -f "$d/expression.rds" ]]; then
          echo "$d"
        fi
      done > "$MANIFEST_FILE"
  N=$(wc -l < "$MANIFEST_FILE" | tr -d ' ')
  if [[ "$N" -eq 0 ]]; then
    echo "ERROR: no bundleable folders under $STUDIES_ABS" >&2
    exit 1
  fi
  echo "Found $N studies under $STUDIES_ABS"
fi
echo "Manifest -> $MANIFEST_FILE"
echo

# Per-task array submission -------------------------------------------------
ARRAY_JOBNAME="scminer_portal_${STAMP}[1-${N}]"
ARRAY_BSUB=(
  bsub
    -J "$ARRAY_JOBNAME"
    -P "$PROJECT"
    -q "$QUEUE"
    -n "$CORES"
    -R "rusage[mem=${MEM}]"
    -W "$WALL"
    -o "paper/logs/portal_studies_%J_%I.out"
    -e "paper/logs/portal_studies_%J_%I.err"
    -env "all, MANIFEST_FILE=${MANIFEST_FILE}"
)

# Merge job depends on the entire array completing successfully ------------
MERGE_BSUB=(
  bsub
    -J "scminer_portal_merge_${STAMP}"
    -P "$PROJECT"
    -q "$QUEUE"
    -n 1
    -R "rusage[mem=4000]"
    -W "0:30"
    -w "done(scminer_portal_${STAMP})"
    -o "paper/logs/portal_merge_${STAMP}.out"
    -e "paper/logs/portal_merge_${STAMP}.err"
)

MERGE_CMD='Rscript paper/portal_merge.R paper/metrics/portal_studies.tsv paper/metrics/portal_studies_*.tsv'

if [[ "$DRY_RUN" == "1" ]]; then
  echo "[dry-run] array submission:"
  printf '    %q ' "${ARRAY_BSUB[@]}"; echo "< paper/portal_studies.bsub"
  echo "[dry-run] merge submission:"
  printf '    %q ' "${MERGE_BSUB[@]}"; echo "\"$MERGE_CMD\""
  exit 0
fi

echo "Submitting array job: ${ARRAY_JOBNAME}"
"${ARRAY_BSUB[@]}" < paper/portal_studies.bsub

echo
echo "Submitting merge job (runs when array completes):"
"${MERGE_BSUB[@]}" "$MERGE_CMD"

echo
echo "Watch progress:"
echo "    bjobs -A         # array status"
echo "    bjobs -J scminer_portal_${STAMP}"
echo "    tail -f paper/logs/portal_studies_*_*.out"
echo
echo "Final merged TSV will appear at paper/metrics/portal_studies.tsv"
