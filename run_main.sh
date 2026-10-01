#!/usr/bin/env bash
# Full pipeline driver for Orthohantavirus reassortment analysis.
# Run from the project root via Git Bash:  bash run_main.sh
#
# Usage:
#   bash run_main.sh              # full run (~22-26 h)
#   bash run_main.sh --preflight  # archive + Step 1 only, then stop to let you inspect
#
# What it does:
#   0. Pre-flight checks (tools, data package, git state, tool versions)
#   1. Archives output/ and data/rpi_calibration.json to a sibling directory
#      OUTSIDE the repo so no pipeline step can accidentally read stale files
#   2. Runs all 5 pipeline stages in order, aborting on any failure
#   3. Writes STATUS.log in the project root with per-step timestamps and commit hash
#
# Archive location: ../  _pipeline_archive/run_YYYYMMDD_HHMMSS/
# Seed 20260904 matches Ortho_rerun08 for rank comparisons.

set -euo pipefail

PREFLIGHT=false
for arg in "$@"; do
    [[ "$arg" == "--preflight" ]] && PREFLIGHT=true
done

# ── Paths ──────────────────────────────────────────────────────────────────
PROJ="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MAFFT="$PROJ/tools/mafft-win/mafft.bat"
IQTREE="$PROJ/tools/iqtree-3.1.3-Windows/bin/iqtree3.exe"
PYTHON="python"
LOG="$PROJ/STATUS.log"

# ── Helpers ────────────────────────────────────────────────────────────────
ts()  { date "+%Y-%m-%dT%H:%M:%S"; }
log() { echo "[$(ts)] $*" | tee -a "$LOG"; }
die() { log "FATAL: $*"; exit 1; }

step() {
    local name="$1"; shift
    log "=== START $name ==="
    "$@" 2>&1 | tee -a "$LOG"
    log "=== DONE  $name ==="
}

# ── Pre-flight checks ──────────────────────────────────────────────────────
cd "$PROJ"
: > "$LOG"   # truncate/create STATUS.log

log "Orthohantavirus pipeline — full scratch run"
log "Project root : $PROJ"
log "Git commit   : $(git rev-parse HEAD)"
log "Git branch   : $(git rev-parse --abbrev-ref HEAD)"
log "Python       : $($PYTHON --version 2>&1)"
[[ "$PREFLIGHT" == true ]] && log "Mode         : PREFLIGHT (archive + Step 1 only)"

# Tools must exist as files
[[ -f "$MAFFT"  ]] || die "MAFFT not found at $MAFFT"
[[ -f "$IQTREE" ]] || die "IQ-TREE not found at $IQTREE"

# Verify tools actually run (not just present)
log "IQ-TREE version: $("$IQTREE" --version 2>&1 | head -1)"
MAFFT_RESPONSE="$("$MAFFT" --help 2>&1 | head -3)"
log "MAFFT response : $(echo "$MAFFT_RESPONSE" | tr '\n' ' ')"
echo "$MAFFT_RESPONSE" | grep -qi "mafft\|environment\|preparing" \
    || die "MAFFT did not produce expected output — check $MAFFT"

# NCBI data package must be present
[[ -d "$PROJ/all_hanta" ]] || die "NCBI data package all_hanta/ not found under $PROJ"

# Warn if working tree is dirty (don't abort — user may have uncommitted README edits etc.)
if [[ -n "$(git status --porcelain)" ]]; then
    log "WARNING: working tree is not clean:"
    git status --short 2>&1 | tee -a "$LOG"
fi

# ── Archive existing outputs — outside the repo ────────────────────────────
# Archive goes to ../  _pipeline_archive/run_YYYYMMDD_HHMMSS/ so no pipeline
# step can accidentally glob or walk into it.
ARCHIVE_ROOT="$PROJ/../_pipeline_archive/run_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$ARCHIVE_ROOT"
log "Archive root : $ARCHIVE_ROOT"

if [[ -d "$PROJ/output" ]]; then
    log "Archiving output/ → $ARCHIVE_ROOT/output"
    mv "$PROJ/output" "$ARCHIVE_ROOT/output"
fi
mkdir -p "$PROJ/output"

# Only rpi_calibration.json is pipeline-produced; country_centroids.csv and
# incidental_hosts.txt are static inputs that must stay in place.
if [[ -f "$PROJ/data/rpi_calibration.json" ]]; then
    log "Archiving data/rpi_calibration.json → $ARCHIVE_ROOT/rpi_calibration.json"
    cp "$PROJ/data/rpi_calibration.json" "$ARCHIVE_ROOT/rpi_calibration.json"
    # Leave the original in place — --recalibrate overwrites it anyway.
    # (cp rather than mv so the run doesn't fail if rpi_core errors before writing a new one)
fi

# ══════════════════════════════════════════════════════════════════════════
# Step 1 — Build isolate pool  (~5 min)
# ══════════════════════════════════════════════════════════════════════════
step "Step1_prepare_dataset" \
    $PYTHON "$PROJ/prepare_dataset.py" \
        --input  "$PROJ" \
        --output "$PROJ/output"

# ── Preflight stop ─────────────────────────────────────────────────────────
if [[ "$PREFLIGHT" == true ]]; then
    log ""
    log "=== PREFLIGHT COMPLETE ==="
    log "Inspect output/ — expected files:"
    log "  output/phylo_rpi_ga_isolates.json"
    log "  output/phylo_rpi_ga_genomic.json"
    log "  output/ga_pool_terminal_complete.csv"
    log "  output/ga_pool_strict_both.csv"
    log "  output/ga_pools_manifest.json"
    ls -lh "$PROJ/output/" 2>&1 | tee -a "$LOG" || true
    log ""
    log "If output looks correct, start the full run:"
    log "  bash run_main.sh"
    log ""
    log "Full log: $LOG"
    exit 0
fi

# ══════════════════════════════════════════════════════════════════════════
# Step 2 — MAFFT + IQ-TREE phylogenies  (~18-22 h)
# ══════════════════════════════════════════════════════════════════════════
step "Step2_phylogeny" \
    $PYTHON "$PROJ/phylogeny.py" \
        --output     "$PROJ/output" \
        --mafft-bin  "$MAFFT" \
        --iqtree-bin "$IQTREE" \
        --phylo-seed 42 \
        --phylo-scope both

# ══════════════════════════════════════════════════════════════════════════
# Step 3a — Recalibrate fallback values, score all triplets  (~1 min)
# Writes fresh data/rpi_calibration.json (overwrites the archived copy's original)
# ══════════════════════════════════════════════════════════════════════════
step "Step3a_rpi_core_run" \
    $PYTHON "$PROJ/rpi_core.py" run \
        --recalibrate \
        --reference-input "$PROJ"

# ══════════════════════════════════════════════════════════════════════════
# Step 3b — GA search  (~2-4 h)
# ══════════════════════════════════════════════════════════════════════════
step "Step3b_ga_search" \
    $PYTHON "$PROJ/Algorithm.py" \
        --reference-input "$PROJ/output" \
        --output          "$PROJ/output/ga_search" \
        --population 150 \
        --generations 100 \
        --seed 20260904 \
        --with-phylogeny

# ══════════════════════════════════════════════════════════════════════════
# Step 4 — Robustness sweep  (~5-10 min)
# ══════════════════════════════════════════════════════════════════════════
step "Step4_robustness_sweep" \
    $PYTHON "$PROJ/rpi_core.py" robustness-sweep \
        --reference-input "$PROJ" \
        --output "$PROJ/output/rpi_parameter_robustness_sweep.json" \
        --figure "$PROJ/output/figures/parameter_robustness_sweep.png"

# ══════════════════════════════════════════════════════════════════════════
# Step 5 — Figures
# ══════════════════════════════════════════════════════════════════════════

step "Step5a_figures_before" \
    $PYTHON "$PROJ/figures.py" before \
        --input  "$PROJ/output" \
        --output "$PROJ/output/figures/before"

# after figures for both pools (skip gracefully if a pool JSON is absent)
for POOL in terminal_complete strict_both; do
    GA_JSON="$PROJ/output/ga_search/$POOL/seed_20260904.json"
    if [[ -f "$GA_JSON" ]]; then
        step "Step5b_figures_after_$POOL" \
            $PYTHON "$PROJ/figures.py" after \
                --input          "$PROJ/output" \
                --ga-search-json "$GA_JSON" \
                --output         "$PROJ/output/figures/after_$POOL"
    else
        log "SKIP Step5b after_$POOL: $GA_JSON not found"
    fi
done

step "Step5c_figures_trees" \
    $PYTHON "$PROJ/figures.py" trees \
        --input  "$PROJ/output" \
        --output "$PROJ/output/figures/trees" \
        --format pdf

for SEG in S M L; do
    step "Step5d_figures_overview_$SEG" \
        $PYTHON "$PROJ/figures.py" overview \
            --input   "$PROJ/output" \
            --output  "$PROJ/output/figures/overview" \
            --segment "$SEG"
done

# ══════════════════════════════════════════════════════════════════════════
# Done
# ══════════════════════════════════════════════════════════════════════════
log ""
log "=== PIPELINE COMPLETE ==="
log ""
log "Post-run diagnostics — run manually after inspecting results:"
log "  # How many of the top-10 winners have neutral factors (fallback-affected)?"
log "  python -c \""
log "  import csv; rows = sorted(csv.DictReader(open('output/ga_search/terminal_complete/all_scored_seed_20260904.csv')), key=lambda r: float(r.get('fitness') or 0), reverse=True)[:10]"
log "  print(sum(1 for r in rows if r.get('rpi_neutral_factors')), 'of 10 top winners have neutral factors')"
log "  \""
log ""
log "  # Compare winner ranks against archived Ortho_rerun08 run (if present):"
log "  # ls $ARCHIVE_ROOT"
log ""
log "Full log: $LOG"
log "Archive : $ARCHIVE_ROOT"
