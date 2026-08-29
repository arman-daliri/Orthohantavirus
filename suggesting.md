### 1

Suggestion:
since your README notes that runtime "depends mostly on IQ-TREE's ModelFinder step and its ultrafast bootstrap step," you'd remove one of the two dominant costs. But there are important nuances before you do it.

## What you'd actually save

ModelFinder tests dozens of candidate models (GTR, SYM, HKY, ... × +I, +Γ, +I+Γ, +R variants), each requiring likelihood optimization over your alignment. For small alignments (4–20 taxa), this is often the *majority* of wall time. For your big species (PUUV, TULV, HTNV with hundreds of isolates), the tree search + 1,000 UFBoot replicates dominate instead — so the speedup shrinks exactly where you'd want it most.

## What you'd risk

- **Model misspecification.** Your within-species alignments are mostly near-identical sequences — for such data, ModelFinder with BIC often picks *simpler* models (HKY+I, K80, even JC-like). Forcing GTR+I+Γ everywhere can over-parameterize, and GTR+I+Γ specifically has a known identifiability problem (the +I and +Γ parameters trade off against each other).
- **Scope asymmetry.** Species-scope trees (shallow) tolerate mild misspecification well; phylogroup-scope trees (deep, cross-species) are where model choice actually bites. A single fixed model is weakest exactly at your interspecies-reassortment layer.
- **Clarification on what "fixed" means:** `-m GTR+G` fixes the model *form* but still re-estimates all parameters per alignment — so you keep per-group rate heterogeneity; you only skip selection. That's the safe kind of shortcut.

## Recommended middle paths (in order of preference)

1. **Restrict, don't fix.** Use `-mset` to limit ModelFinder to a small relevant set (e.g., `-mset GTR,HKY,SYM,TIM` — with rate variants). Cuts ModelFinder cost substantially while keeping per-alignment selection. Smallest rigor loss.
2. **Select once per segment, then fix.** Run ModelFinder on the largest/most diverse alignment per segment (e.g., PUUV L), adopt its winner (commonly something like `GTR+F+I+G4`), fix it for the remaining groups of that segment — and **document it plus a sensitivity check**: re-run 2–3 species with `-m MFP` and compare RF distances and your reassortment calls. If calls don't change, the shortcut is empirically justified in your own data.
3. **Fix globally without checking.** Fastest, but I'd advise against it for a paper — it's the one version a reviewer can legitimately challenge.

## One thing not to touch

Keep `-B 1000`. Your entire reassortment gate rests on the ≥95 UFBoot threshold, and UFBoot is only calibrated at ≥1,000 replicates. If runtime is still painful after the model shortcut, the honest lever is `--phylo-threads` (already cores−1) or subsampling large species — not fewer bootstraps.

Bottom line: option 1 or 2 gives you most of the speedup with a defensible Methods sentence ("model selection restricted to a candidate set / model fixed per segment after ModelFinder evaluation, with topology and reassortment calls verified unchanged under full `-m MFP` for a subset of species").