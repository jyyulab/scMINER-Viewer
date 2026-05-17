#!/usr/bin/env bash
# paper/benchmarks/run_figure1.sh
#
# Runs the synthetic-sweep benchmark that produces figure 1 (and
# table 2). Sequentially executes:
#
#   1.  Rscript paper/benchmarks/figures.R     7x4 = 28 configs
#   2.  Rscript paper/benchmarks/tables.R      regenerate table 2
#
# Two modes:
#
#   * default     run locally (laptop / HPC interactive node). Logs are
#                 tee'd to paper/logs/figure1_<timestamp>.log.
#   * --hpc       submit as a single LSF bsub job. Useful when the
#                 10K x 10K configs would OOM your laptop.
#
# Usage:
#   ./paper/benchmarks/run_figure1.sh
#   ./paper/benchmarks/run_figure1.sh --hpc --mem 32000 --wall 4:00
#   ./paper/benchmarks/run_figure1.sh --hpc --queue priority --mem 64000
#
# Flags (only meaningful with --hpc, ignored locally):
#   --hpc                 Submit via bsub instead of running locally.
#   --mem <MB>            Memory request (default: 32000 = 32 GB).
#   --cores <n>           Cores per job (default: 4).
#   --wall <hh:mm>        Wall-clock limit (default: 4:00).
#   --queue <name>        LSF queue (default: standard).
#   --project <id>        LSF charge code (default: scminer).
#   --dry-run             Print the bsub command without submitting.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

HPC=0
MEM=32000
CORES=4
WALL="4:00"
QUEUE="standard"
PROJECT="scminer"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --hpc)        HPC=1; shift ;;
    --mem)        MEM="$2"; shift 2 ;;
    --cores)      CORES="$2"; shift 2 ;;
    --wall|-W)    WALL="$2"; shift 2 ;;
    --queue)      QUEUE="$2"; shift 2 ;;
    --project)    PROJECT="$2"; shift 2 ;;
    --dry-run)    DRY_RUN=1; shift ;;
    -h|--help)
      sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "unknown flag: $1" >&2; exit 2 ;;
  esac
done

mkdir -p paper/figures paper/metrics paper/tables paper/logs

STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="paper/logs/figure1_${STAMP}.log"

run_local() {
  echo "scMINER Viewer -- figure 1 synthetic-sweep benchmark"
  echo "  root:    $ROOT"
  echo "  log:     $LOG"
  echo "  start:   $(date)"
  echo

  {
    echo "=== figures.R ==="
    Rscript paper/benchmarks/figures.R
    echo
    echo "=== tables.R ==="
    Rscript paper/benchmarks/tables.R
  } 2>&1 | tee "$LOG"

  echo
  echo "Done."
  echo "  figures: paper/figures/figure2_{A..F}_*.{pdf,png}"
  echo "  table:   paper/tables/figure2_scaling.{md,tsv}"
  echo "  metrics: paper/metrics/{bundle_scaling,discover_scaling,real_study}.tsv"
  echo "  log:     $LOG"
}

submit_hpc() {
  # Inline script so we can `module load R` and chain figures.R + tables.R.
  local inline
  inline="$(cat <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "\${LS_SUBCWD:-$ROOT}"
module load R/4.2.2-rhel8 2>/dev/null || true
module load hdf5 2>/dev/null || module load hdf5/1.10.7 2>/dev/null || module load hdf5/1.10.8 2>/dev/null || true
echo "[bsub] figure1 job=\${LSB_JOBID:-?} host=\$(hostname)"
echo "[bsub] start: \$(date)"
echo
echo "=== figures.R ==="
Rscript paper/benchmarks/figures.R
echo
echo "=== tables.R ==="
Rscript paper/benchmarks/tables.R
echo
echo "[bsub] done -- artifacts under paper/figures/, paper/tables/, paper/metrics/"
EOF
)"

  local cmd=(
    bsub
      -J "scminer_figure1_${STAMP}"
      -P "$PROJECT"
      -q "$QUEUE"
      -n "$CORES"
      -R "rusage[mem=${MEM}]"
      -M "${MEM}"
      -W "$WALL"
      -o "paper/logs/figure1_${STAMP}.out"
      -e "paper/logs/figure1_${STAMP}.err"
  )

  echo "Submitting figure 1 benchmark to LSF"
  echo "  queue:    $QUEUE"
  echo "  cores:    $CORES"
  echo "  mem (MB): $MEM"
  echo "  wall:     $WALL"
  echo

  if [[ "$DRY_RUN" == "1" ]]; then
    echo "[dry-run]:"
    printf '    %q ' "${cmd[@]}"; echo " <<inline-script>>"
    echo
    echo "[dry-run] inline script:"
    printf '%s\n' "$inline" | sed 's/^/    /'
    exit 0
  fi

  printf '%s\n' "$inline" | "${cmd[@]}"

  echo
  echo "Watch progress:"
  echo "    bjobs -J scminer_figure1_${STAMP}"
  echo "    tail -f paper/logs/figure1_${STAMP}.out"
  echo
  echo "When it finishes:"
  echo "    ls   paper/figures/figure2_*"
  echo "    cat  paper/tables/figure2_scaling.md"
}

if [[ "$HPC" == "1" ]]; then
  submit_hpc
else
  run_local
fi
