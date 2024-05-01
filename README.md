# Basophil-immune-model
Code repository for immune models to study the effect of basophils on enhancing cancer cell killing by CD8 T cells.

run the "driver_May2024.m" file, which calls individual model files "core_<filename>.m", calculates percent cancer cell death, and plots results.

core_base.m - contains ODEs for base model of cancer cells, Tregs, and CD8 T cells (no basophils)

core_TregDeath.m - contains ODEs for model where basophils promote Treg death

core_CD8supp.m - contains ODEs for model where basophils promote Treg-mediated suppression of CD8 T cells

core_TregDeath_CD8supp.m - contains ODEs for model where basophils promote death of Tregs and Treg-mediated suppression of CD8 T cells


Will update upon revision and publication.
