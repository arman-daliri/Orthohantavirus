#!/usr/bin/env bash
# Full pipeline driver for Orthohantavirus reassortment analysis.
# Run from the project root:  bash run_main.sh
#
# What it does:
#   1. Archives the current output/ directory (if any) so this run starts clean
#   2. Runs all 5 pipeline steps in order, aborting on any failure
#   3. Writes STATUS.log in the project root with per-step timestamps
#
# Expected wall time: ~20-26 hours (Step 2 dominates)
# Seed 20260904 matches Ortho_rerun08 for rank comparisons.

set -euo pipefail

# ── Paths ──────────────────────────────────────────────────────────────────
PROJ="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MAFFT="$PROJ/tools/mafft-win/mafft.bat"
IQTREE="$PROJ/tools/iqtree-3.1.3-Windows/bin/iqtree3.exe"
PYTHON="python"
LOG="$PROJ/STATUS.log"

# ── Helpers ────────────────────────────────────────────────────────────────
ts() { date "+%Y-%m-%dT%H:%M:%S"; }
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

[[ -f "$MAFFT"  ]] || die "MAFFT not found at $MAFFT"
[[ -f "$IQTREE" ]] || die "IQ-TREE not found at $IQTREE"
[[ -d "$PROJ/all_hanta" ]] || die "NCBI data package all_hanta/ not found"

if [[ -n "$(git status --porcelain)" ]]; then
    log "WARNING: working tree is not clean — uncommitted changes present"
    log "$(git status --short)"
fi

# ── Archive existing output/ ───────────────────────────────────────────────
if [[ -d "$PROJ/output" ]]; then
    ARCHIVE="$PROJ/output/_archive_$(date +%Y%m%d_%H%M%S)"
    log "Archiving existing output/ → $ARCHIVE"
    mv "$PROJ/output" "$ARCHIVE"
fi
mkdir -p "$PROJ/output"

# ══════════════════════════════════════════════════════════════════════════
# Step 1 — Build isolate pool  (~5 min)
# ══════════════════════════════════════════════════════════════════════════
step "Step1_prepare_dataset" \
    $PYTHON "$PROJ/prepare_dataset.py" \
        --input "$PROJ" \
        --output "$PROJ/output"

# ══════════════════════════════════════════════════════════════════════════
# Step 2 — MAFFT + IQ-TREE phylogenies  (~18-22 h)
# ══════════════════════════════════════════════════════════════════════════
step "Step2_phylogeny" \
    $PYTHON "$PROJ/phylogeny.py" \
        --output "$PROJ/output" \
        --mafft-bin  "$MAFFT" \
        --iqtree-bin "$IQTREE" \
        --phylo-seed 42 \
        --phylo-scope both

# ══════════════════════════════════════════════════════════════════════════
# Step 3a — Recalibrate + score all triplets  (~1 min)
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
        --output         "$PROJ/output/ga_search" \
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
log "=== PIPELINE COMPLETE ==="
log ""
log "Post-run diagnostics (run manually after inspecting results):"
log "  Top-10 winners in terminal_complete with non-empty rpi_neutral_factors:"
log "    python -c \""
log "    import json, csv"
log "    rows = list(csv.DictReader(open('output/ga_search/terminal_complete/all_scored_seed_20260904.csv')))"
log "    top = sorted(rows, key=lambda r: float(r.get('fitness') or 0), reverse=True)[:10]"
log "    print(sum(1 for r in top if r.get('rpi_neutral_factors')),'of 10 top winners have neutral factors')"
log "    \""
log ""
log "Full log: $LOG"
