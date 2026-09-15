# Basophil-immune-model
Code repository for computational models to study the effect of basophils on enhancing cancer cell killing by CD8 T cells.

run the "driver.m" file, which calls individual model files "core_FILENAME.m", calculates percent cancer cell death, and plots the simulated species' dynamics.

core_base.m - contains ODEs for base model of cancer cells, Tregs, and CD8 T cells (no basophils)

core_TregDeath.m - contains ODEs for model where basophils promote Treg death

core_CD8supp.m - contains ODEs for model where basophils inhibit Treg-mediated suppression of CD8 T cells

core_TregDeath_CD8prolif.m - contains ODEs for model where basophils promote death of Tregs and inhibits Treg-mediated suppression of CD8 T cells
