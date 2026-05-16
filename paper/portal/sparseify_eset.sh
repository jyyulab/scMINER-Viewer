#!/usr/bin/env bash
# paper/portal/sparseify_eset.sh
#
# One-time bsub wrapper around paper/portal/sparseify_eset.R: converts a
# dense-backed Biobase::ExpressionSet rds into a dgCMatrix-backed one,
# in a single high-mem LSF job. Use this to "prime" any study whose
# raw expression.rds blows past the standard 200 GB ceiling during
# readRDS (Covid650k / 2317 is the canonical case).
#
# Usage:
#   ./paper/portal/sparseify_eset.sh \
#       --in  /research/.../Covid650k/expression.rds
#
#   ./paper/portal/sparseify_eset.sh \
#       --in  /research/.../Covid650k/expression.rds \
#       --out /research/.../Covid650k/expression.sparse.rds \
#       --mem 400000 --wall 4:00
#
#   # Verify only -- no write:
#   ./paper/portal/sparseify_eset.sh \
#       --in  /research/.../Covid650k/expression.rds \
#       --verify-only
#
# Flags:
#   --in  <path>          (required)  source .rds (ExpressionSet)
#   --out <path>          destination .rds (default: <stem>.sparse.rds)
#   --force               overwrite existing destination
#   --verify-only         load + report only; do not write
#   --mem  <MB>           memory request (default: 400000 = 400 GB)
#   --cores <n>           cores (default: 1)
#   --wall  <hh:mm>       wall-clock (default: 6:00)
#   --queue <name>        LSF queue (default: standard)
#   --project <id>        LSF charge code (default: scminer)
#   --dry-run             print the bsub command without submitting

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

SRC=""
DST=""
FORCE=0
VERIFY_ONLY=0
MEM=400000
CORES=1
WALL="6:00"
QUEUE="standard"
PROJECT="scminer"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --in)            SRC="$2"; shift 2 ;;
    --out)           DST="$2"; shift 2 ;;
    --force)         FORCE=1; shift ;;
    --verify-only)   VERIFY_ONLY=1; shift ;;
    --mem)           MEM="$2"; shift 2 ;;
    --cores)         CORES="$2"; shift 2 ;;
    --wall|-W)       WALL="$2"; shift 2 ;;
    --queue)         QUEUE="$2"; shift 2 ;;
    --project)       PROJECT="$2"; shift 2 ;;
    --dry-run)       DRY_RUN=1; shift ;;
    -h|--help)
      sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "unknown flag: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$SRC" ]]; then
  echo "ERROR: --in <path> is required" >&2
  exit 2
fi
if [[ ! -f "$SRC" ]]; then
  echo "ERROR: source does not exist: $SRC" >&2
  exit 2
fi

mkdir -p paper/logs

# Build the inner Rscript command --------------------------------------------
R_ARGS=(--in "$SRC")
if [[ -n "$DST" ]];           then R_ARGS+=(--out "$DST"); fi
if [[ "$FORCE" == "1" ]];     then R_ARGS+=(--force); fi
if [[ "$VERIFY_ONLY" == "1" ]]; then R_ARGS+=(--verify-only); fi

STAMP="$(date +%Y%m%d_%H%M%S)"
SRC_TAG="$(basename "${SRC%.*}")"

# Use a small inline bsub script so we can `module load R/4.2.2-rhel8`
# and run the Rscript with the right env on whichever compute node.
INLINE_SCRIPT="$(cat <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "\${LS_SUBCWD:-$ROOT}"
module load R/4.2.2-rhel8 2>/dev/null || true
module load hdf5 2>/dev/null || module load hdf5/1.10.7 2>/dev/null || module load hdf5/1.10.8 2>/dev/null || true
echo "[bsub] sparseify_eset job=\${LSB_JOBID:-?} host=\$(hostname)"
echo "[bsub] src: $SRC"
echo "[bsub] dst: ${DST:-<auto>}"
Rscript paper/portal/sparseify_eset.R ${R_ARGS[*]@Q}
EOF
)"

BSUB_CMD=(
  bsub
    -J "sparseify_${SRC_TAG}_${STAMP}"
    -P "$PROJECT"
    -q "$QUEUE"
    -n "$CORES"
    -R "rusage[mem=${MEM}]"
    -M "${MEM}"
    -W "$WALL"
    -o "paper/logs/sparseify_${SRC_TAG}_${STAMP}.out"
    -e "paper/logs/sparseify_${SRC_TAG}_${STAMP}.err"
)

echo "Submitting sparseify job"
echo "  source:    $SRC"
echo "  dest:      ${DST:-<auto: <stem>.sparse.rds>}"
echo "  queue:     $QUEUE"
echo "  cores:     $CORES"
echo "  mem (MB):  $MEM"
echo "  wall:      $WALL"
echo "  R args:    ${R_ARGS[*]}"
echo

if [[ "$DRY_RUN" == "1" ]]; then
  echo "[dry-run]"
  printf '    %q ' "${BSUB_CMD[@]}"; echo " <<inline-script>>"
  echo
  echo "[dry-run] inline script:"
  printf '%s\n' "$INLINE_SCRIPT" | sed 's/^/    /'
  exit 0
fi

printf '%s\n' "$INLINE_SCRIPT" | "${BSUB_CMD[@]}"

echo
echo "Watch progress:"
echo "    bjobs -J sparseify_${SRC_TAG}_${STAMP}"
echo "    tail -f paper/logs/sparseify_${SRC_TAG}_${STAMP}.out"
echo
echo "When it finishes, point the YAML at the new file:"
echo "    yq -y '.input.expression = \"${DST:-<auto>.sparse.rds>}\"' \\"
echo "        -i paper/configs/<id>.yaml"
