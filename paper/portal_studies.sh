#!/usr/bin/env bash
# paper/portal_studies.sh
#
# Convenience wrapper for paper/portal_studies.R. Use this when running
# directly (no scheduler) on an HPC interactive node or local machine.
# For batch submission on LSF, see paper/portal_studies.bsub.
#
# Usage:
#   ./paper/portal_studies.sh --studies-root /path/to/studies
#   ./paper/portal_studies.sh --studies-root <path> --cap 10000
#   ./paper/portal_studies.sh --studies-root <path> --only Tex,Tregs
#
# Each study folder under --studies-root must contain at minimum
# expression.rds (ExpressionSet). activity.rds and network.txt are
# optional but recommended for complete metrics.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

mkdir -p paper/metrics paper/logs

LOG="paper/logs/portal_studies_$(date +%Y%m%d_%H%M%S).log"

echo "scMINER Viewer -- portal-study benchmark (real data)"
echo "  root:    $ROOT"
echo "  log:     $LOG"
echo "  args:    $*"
echo

Rscript paper/portal_studies.R "$@" 2>&1 | tee "$LOG"

echo
echo "Done. Metrics: paper/metrics/portal_studies.tsv"
echo "Log:           $LOG"
