# Basophil-immune-model
Code repository for immune models to study the effect of basophils on enhancing cancer cell killing by CD8 T cells.

run the "driver_NEW.m" file, which calls individual model files "core_<filename>.m", calculates percent cancer cell death, and plots results.

core_base.m - contains ODEs for base model of cancer cells, Tregs, and CD8 T cells (no basophils)
core_TregDeath.m - contains ODEs for model where basophils promote Treg death
core_CD8prolif.m - contains ODEs for model where basophils promote CD8 T cell proliferation
core_TregDeath_CD8prolif.m - contains ODEs for model where basophils promote death of Tregs and CD8 T cell proliferation

Will update upon revision and publication.
