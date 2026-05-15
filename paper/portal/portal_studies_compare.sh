#!/usr/bin/env bash
# paper/portal/portal_studies_compare.sh
#
# With-vs-without-TFsig benchmark fan-out for the 29 portal studies.
#
# Submits TWO bsub job arrays back-to-back:
#
#   1. expression-only  (no activity.rds, no networks.txt)
#   2. full             (expression + activity + networks; current default)
#
# The 'full' array depends on `ended(expression-only)` so it starts only
# after every expression-only task has finished -- this keeps the two
# modes from racing on the same compute-node page cache and the same
# parallel-fs read bandwidth, which would smear the wall-time and
# warm-load measurements.
#
# A third "compare" job depends on `ended(full)` and produces
#   paper/metrics/portal_studies_compare.tsv
# by joining the two per-mode TSVs on studyID.
#
# Each array task reuses paper/portal/portal_studies.bsub with MODE set,
# so per-study TSVs land at:
#   paper/metrics/portal_studies_<id>_expression-only.tsv
#   paper/metrics/portal_studies_<id>_full.tsv
#
# Usage (from project root):
#   ./paper/portal/portal_studies_compare.sh --configs-dir paper/configs
#   ./paper/portal/portal_studies_compare.sh --configs-dir paper/configs \
#       --mem 64000 --wall 8:00 --queue "standard priority"
#   ./paper/portal/portal_studies_compare.sh --configs-dir paper/configs \
#       --only 2317,2327 --dry-run
#
# Flags:
#   --configs-dir <dir>    Folder of YAML configs (default: paper/configs).
#   --queue "<q1 q2 ...>"  LSF queue list (default: "standard priority").
#                          Tasks within a mode are distributed across queues.
#   --mem <MB>             Memory request per task (default: 64000).
#   --cores <n>            Cores per task (default: 4).
#   --wall <hh:mm>         Wall-clock limit per task (default: 8:00).
#   --project <id>         Charge code (default: scminer).
#   --only <csv>           Comma-separated study IDs to limit to (e.g. 2317,2327).
#                          Useful for re-running a subset after fixing data.
#   --dry-run              Print the bsub commands without submitting.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

CONFIGS_DIR="paper/configs"
QUEUE="standard priority"
MEM=64000
CORES=4
WALL="8:00"
PROJECT="scminer"
ONLY=""
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --configs-dir)  CONFIGS_DIR="$2"; shift 2 ;;
    --queue)        QUEUE="$2"; shift 2 ;;
    --mem)          MEM="$2"; shift 2 ;;
    --cores)        CORES="$2"; shift 2 ;;
    --wall|-W)      WALL="$2"; shift 2 ;;
    --project)      PROJECT="$2"; shift 2 ;;
    --only)         ONLY="$2"; shift 2 ;;
    --dry-run)      DRY_RUN=1; shift ;;
    -h|--help)
      sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "unknown flag: $1" >&2; exit 2 ;;
  esac
done

if [[ ! -d "$CONFIGS_DIR" ]]; then
  echo "ERROR: --configs-dir does not exist: $CONFIGS_DIR" >&2
  exit 2
fi

mkdir -p paper/metrics paper/logs paper/hpc

STAMP="$(date +%Y%m%d_%H%M%S)"
MANIFEST_FILE="paper/hpc/manifest_compare_${STAMP}.txt"

# Build the manifest: one absolute YAML path per line. Apply --only by
# filtering on filename stem (2317, 2327, ...) so the array index range
# matches the actual task count.
CONFIGS_ABS="$(cd "$CONFIGS_DIR" && pwd)"
ONLY_SET=""
if [[ -n "$ONLY" ]]; then
  ONLY_SET="$(echo "$ONLY" | tr ',' '|')"
fi
if [[ -n "$ONLY_SET" ]]; then
  find "$CONFIGS_ABS" -mindepth 1 -maxdepth 1 -type f \
       \( -name '*.yaml' -o -name '*.yml' \) \
    | sort \
    | awk -v re="^($ONLY_SET)\\.ya?ml$" '
        { n = split($0, parts, "/"); if (parts[n] ~ re) print }
      ' > "$MANIFEST_FILE"
else
  find "$CONFIGS_ABS" -mindepth 1 -maxdepth 1 -type f \
       \( -name '*.yaml' -o -name '*.yml' \) \
    | sort > "$MANIFEST_FILE"
fi

N=$(wc -l < "$MANIFEST_FILE" | tr -d ' ')
if [[ "$N" -eq 0 ]]; then
  echo "ERROR: no YAML configs matched in $CONFIGS_ABS" \
       "${ONLY:+(filter: $ONLY)}" >&2
  exit 1
fi

echo "Compare benchmark: with-TFsig vs expression-only"
echo "  configs:   $CONFIGS_ABS"
echo "  manifest:  $MANIFEST_FILE  ($N studies)"
echo "  queue(s):  $QUEUE"
echo "  mem (MB):  $MEM   cores: $CORES   wall: $WALL"
echo

# Split tasks across queues per mode -----------------------------------------
# Same trick as portal_studies_hpc.sh: LSF's `-q "q1 q2"` does not
# distribute, so we issue one sub-array per queue with a disjoint index
# range over the same manifest.
read -r -a QUEUE_ARR <<< "$QUEUE"
N_Q="${#QUEUE_ARR[@]}"
if [[ "$N_Q" -lt 1 ]]; then
  echo "ERROR: --queue is empty" >&2
  exit 2
fi
CHUNK=$(( (N + N_Q - 1) / N_Q ))

submit_mode() {
  local mode="$1" depends_on="${2:-}"
  local -a mode_arrays=()
  local q_i start end qname jobname
  echo "==== mode: $mode ===="
  for (( q_i = 0; q_i < N_Q; q_i++ )); do
    start=$(( q_i * CHUNK + 1 ))
    end=$(( (q_i + 1) * CHUNK ))
    if [[ "$end" -gt "$N" ]]; then end="$N"; fi
    if [[ "$start" -gt "$N" ]]; then break; fi
    qname="${QUEUE_ARR[$q_i]}"
    jobname="scminer_cmp_${STAMP}_${mode}_${qname}[${start}-${end}]"
    mode_arrays+=("scminer_cmp_${STAMP}_${mode}_${qname}")
    local cmd=(
      bsub
        -J "$jobname"
        -P "$PROJECT"
        -q "$qname"
        -n "$CORES"
        -R "rusage[mem=${MEM}]"
        -W "$WALL"
        -o "paper/logs/portal_studies_${mode}_%J_%I.out"
        -e "paper/logs/portal_studies_${mode}_%J_%I.err"
        -env "all, MANIFEST_FILE=${MANIFEST_FILE}, MODE=${mode}"
    )
    if [[ -n "$depends_on" ]]; then
      cmd+=( -w "$depends_on" )
    fi
    if [[ "$DRY_RUN" == "1" ]]; then
      printf '  [dry-run] '
      printf '%q ' "${cmd[@]}"
      echo "< paper/portal/portal_studies.bsub"
    else
      echo "  Submitting $jobname on $qname"
      "${cmd[@]}" < paper/portal/portal_studies.bsub
    fi
  done
  # Build "ended(arr1) && ended(arr2) && ..." for the next stage's -w.
  local dep=""
  for nm in "${mode_arrays[@]}"; do
    if [[ -n "$dep" ]]; then dep+=" && "; fi
    dep+="ended(${nm})"
  done
  # Echo the dependency string so the caller can capture it.
  echo "$dep"
}

# Capture stdout from submit_mode but still let its progress lines show.
# We grep the last line, which is the dependency string.
EXPR_OUT="$(submit_mode "expression-only" "" | tee /dev/tty)"
EXPR_DEP="$(printf '%s\n' "$EXPR_OUT" | tail -n 1)"

FULL_OUT="$(submit_mode "full" "$EXPR_DEP" | tee /dev/tty)"
FULL_DEP="$(printf '%s\n' "$FULL_OUT" | tail -n 1)"

# Compare job ----------------------------------------------------------------
# Joins paper/metrics/portal_studies_*_expression-only.tsv and
# paper/metrics/portal_studies_*_full.tsv on studyID, emits
# paper/metrics/portal_studies_compare.tsv.
COMPARE_CMD='Rscript paper/portal/portal_compare.R'
COMPARE_CMD+=' --expr-only-glob "paper/metrics/portal_studies_*_expression-only.tsv"'
COMPARE_CMD+=' --full-glob      "paper/metrics/portal_studies_*_full.tsv"'
COMPARE_CMD+=' --out            "paper/metrics/portal_studies_compare.tsv"'

COMPARE_BSUB=(
  bsub
    -J "scminer_cmp_${STAMP}_compare"
    -P "$PROJECT"
    -q "${QUEUE_ARR[0]}"
    -n 1
    -R "rusage[mem=4000]"
    -W "0:30"
    -w "$FULL_DEP"
    -o "paper/logs/portal_compare_${STAMP}.out"
    -e "paper/logs/portal_compare_${STAMP}.err"
)

echo
echo "==== compare job ===="
if [[ "$DRY_RUN" == "1" ]]; then
  printf '  [dry-run] '
  printf '%q ' "${COMPARE_BSUB[@]}"
  echo "\"$COMPARE_CMD\""
  exit 0
fi
echo "  Submitting compare job (waits for: $FULL_DEP)"
"${COMPARE_BSUB[@]}" "$COMPARE_CMD"

echo
echo "Watch progress:"
echo "    bjobs -A                                                 # all arrays"
echo "    bjobs -J 'scminer_cmp_${STAMP}_expression-only_*'        # mode 1"
echo "    bjobs -J 'scminer_cmp_${STAMP}_full_*'                   # mode 2"
echo "    bjobs -J 'scminer_cmp_${STAMP}_compare'                  # join"
echo "    tail -f paper/logs/portal_studies_expression-only_*_*.out"
echo "    tail -f paper/logs/portal_studies_full_*_*.out"
echo
echo "Final comparison TSV: paper/metrics/portal_studies_compare.tsv"
