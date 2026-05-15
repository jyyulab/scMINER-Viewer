#!/usr/bin/env bash
# paper/portal/portal_studies_hpc.sh
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
#   ./paper/portal/portal_studies_hpc.sh --configs-dir paper/configs
#   ./paper/portal/portal_studies_hpc.sh --studies-root /research/.../studies
#   ./paper/portal/portal_studies_hpc.sh --configs-dir <path> --mem 64000 -W 12:00
#   ./paper/portal/portal_studies_hpc.sh --configs-dir <path> --queue "standard priority compbio"
#
# Flags:
#   --configs-dir <dir>    Folder of YAML configs (one per study).
#   --studies-root <dir>   Folder of per-study subfolders (legacy mode).
#   --queue "<q1 q2 ...>"  Space-separated LSF queue list; LSF picks
#                          the first with capacity. Quote for >1 queue.
#                          (default: "standard priority")
#   --mem <MB>             Memory request per task (default: 32000).
#   --cores <n>            Cores per task (default: 4).
#   --wall <hh:mm>         Wall-clock limit per task (default: 6:00).
#   --project <id>         Charge code (default: scminer).
#   --dry-run              Print bsub commands without submitting.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

CONFIGS_DIR=""
STUDIES_ROOT=""
# Space-separated list of LSF queues. LSF dispatches each task to the
# first queue in the list with available capacity, so listing multiple
# queues speeds startup when one is congested.
QUEUE="standard priority"
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

# Split tasks across queues -------------------------------------------------
#
# LSF's native multi-queue syntax (`-q "q1 q2"`) does NOT distribute jobs
# across queues — at dispatch time LSF just picks whichever queue has
# capacity, so in practice everything tends to land in the first one.
# To actually spread the array, we submit one sub-array per queue with a
# disjoint index range over the same manifest. The merge job depends on
# every sub-array completing.
read -r -a QUEUE_ARR <<< "$QUEUE"
N_Q="${#QUEUE_ARR[@]}"
if [[ "$N_Q" -lt 1 ]]; then
  echo "ERROR: --queue is empty" >&2
  exit 2
fi
# Chunk size, ceiling division so every task is covered.
CHUNK=$(( (N + N_Q - 1) / N_Q ))

echo "Distributing $N tasks across $N_Q queue(s) (chunk size $CHUNK):"
ARRAY_NAMES=()
ARRAY_PLAN=()  # for dry-run printing: "queue start end"
for (( q_i = 0; q_i < N_Q; q_i++ )); do
  start=$(( q_i * CHUNK + 1 ))
  end=$(( (q_i + 1) * CHUNK ))
  if [[ "$end" -gt "$N" ]]; then end="$N"; fi
  if [[ "$start" -gt "$N" ]]; then break; fi
  qname="${QUEUE_ARR[$q_i]}"
  jobname="scminer_portal_${STAMP}_${qname}[${start}-${end}]"
  ARRAY_NAMES+=("scminer_portal_${STAMP}_${qname}")
  ARRAY_PLAN+=("$qname $start $end")
  printf '  %-12s indices [%d-%d]  (%d task%s)\n' \
      "$qname" "$start" "$end" $((end - start + 1)) \
      "$(if [[ $((end - start + 1)) -ne 1 ]]; then echo s; fi)"
done
echo

# Build the merge job's dependency string: ended(arr1) && ended(arr2) && ...
# We use ended() rather than done() so the merge fires even when some
# tasks fail -- e.g. a study has an unfixable data issue -- and
# concatenates whatever per-study TSVs landed on disk.
MERGE_DEP=""
for nm in "${ARRAY_NAMES[@]}"; do
  if [[ -n "$MERGE_DEP" ]]; then MERGE_DEP+=" && "; fi
  MERGE_DEP+="ended(${nm})"
done

MERGE_CMD='Rscript paper/portal/portal_merge.R paper/metrics/portal_studies.tsv paper/metrics/portal_studies_*.tsv'

submit_array_for_queue() {
  local qname="$1" start="$2" end="$3"
  local jobname="scminer_portal_${STAMP}_${qname}[${start}-${end}]"
  local cmd=(
    bsub
      -J "$jobname"
      -P "$PROJECT"
      -q "$qname"
      -n "$CORES"
      -R "rusage[mem=${MEM}]"
      -W "$WALL"
      -o "paper/logs/portal_studies_%J_%I.out"
      -e "paper/logs/portal_studies_%J_%I.err"
      -env "all, MANIFEST_FILE=${MANIFEST_FILE}"
  )
  if [[ "$DRY_RUN" == "1" ]]; then
    echo "[dry-run] ${qname} sub-array:"
    printf '    %q ' "${cmd[@]}"; echo "< paper/portal/portal_studies.bsub"
  else
    echo "Submitting ${jobname} on queue ${qname}"
    "${cmd[@]}" < paper/portal/portal_studies.bsub
  fi
}

submit_merge_job() {
  local cmd=(
    bsub
      -J "scminer_portal_merge_${STAMP}"
      -P "$PROJECT"
      -q "${QUEUE_ARR[0]}"
      -n 1
      -R "rusage[mem=4000]"
      -W "0:30"
      -w "$MERGE_DEP"
      -o "paper/logs/portal_merge_${STAMP}.out"
      -e "paper/logs/portal_merge_${STAMP}.err"
  )
  if [[ "$DRY_RUN" == "1" ]]; then
    echo "[dry-run] merge submission:"
    printf '    %q ' "${cmd[@]}"; echo "\"$MERGE_CMD\""
  else
    echo "Submitting merge job (waits for: $MERGE_DEP)"
    "${cmd[@]}" "$MERGE_CMD"
  fi
}

for plan in "${ARRAY_PLAN[@]}"; do
  read -r qname start end <<< "$plan"
  submit_array_for_queue "$qname" "$start" "$end"
done
echo
submit_merge_job

if [[ "$DRY_RUN" == "1" ]]; then
  exit 0
fi

echo
echo "Watch progress:"
echo "    bjobs -A                                       # all arrays"
for nm in "${ARRAY_NAMES[@]}"; do
  echo "    bjobs -J '${nm}*'                            # ${nm}"
done
echo "    tail -f paper/logs/portal_studies_*_*.out"
echo
echo "Final merged TSV will appear at paper/metrics/portal_studies.tsv"
