# Longitudinal UAV remote sensing for wheat breeding decisions

Reproducibility materials for the manuscript:

**Longitudinal UAV remote sensing for wheat breeding decisions: cross-experiment transferability and early selection utility**

Target journal: *Remote Sensing Applications: Society and Environment* (RSASE).

## Contents

- `data/CEA_wheat_phenomics_154_genotypes.csv` — genotype-level analytical dataset (154 wheat genotypes; DH, GY, and 120 longitudinal spectral predictors from 10 UAV flights).
- `scripts/RSASE_master_pipeline.py` — production reference for the audited RSASE revision.
- `results/cross_experiment_spearman.csv` — Spearman rank correlations calculated from the final held-out-domain predictions.
- `results/longitudinal_xgboost_fixed_metrics.csv` — longitudinal XGBoost results after freezing the model specification for the RSASE revision.
- `results/earliest_useful_xgboost_fixed.csv` — earliest-useful-window results for the frozen XGBoost specification.

## Experimental domains

- E1: 49 cultivars
- E2: 56 breeding populations
- E3: 19 breeding populations
- E4: 30 elite lines

All four experiments were conducted in the same site-year and are treated as contrasting experimental domains, not as external multi-environment validation.

## Frozen longitudinal XGBoost specification

The RSASE revision uses one fixed XGBoost configuration across cumulative windows F1-F10:

- n_estimators = 120
- max_depth = 2
- learning_rate = 0.05
- subsample = 0.80
- colsample_bytree = 0.80
- min_child_weight = 1
- reg_lambda = 1.0
- reg_alpha = 0.0
- gamma = 0.0
- objective = reg:squarederror
- random_state = 20260922
- n_jobs = 1

This fixed specification is used to isolate temporal information accumulation rather than re-optimizing model complexity at every flight window.

## Cross-experiment metrics

Cross-experiment reporting now separates:

- Pearson r — linear predictive concordance;
- Spearman rho — ranking transferability;
- R2 — absolute calibration;
- RMSE and MAE — prediction error.

## Prior related publication

The E1 field campaign overlaps with the 49-genotype longitudinal dataset previously analyzed in:

Silva et al. (2026). *Optimizing indirect selection of tropical wheat genotypes using high-throughput longitudinal phenotyping and trait relationships*. Remote Sensing Applications: Society and Environment, 41, 101892. https://doi.org/10.1016/j.rsase.2026.101892

That paper addressed longitudinal repeatability, trait relationships, and indirect selection. The present study addresses distinct questions: leakage-safe prediction, cross-experiment transferability, cumulative information timing, and stability of selection decisions.

## Reproducibility note

The Python workflow is the production reference for the RSASE revision. The previous CEA folder is retained only as an archive of the earlier submission-stage workflow.

Corresponding author: **Maicon Nardino**, Federal University of Viçosa, Brazil.
