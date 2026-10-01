# Longitudinal UAV remote sensing for wheat breeding decisions

Reproducibility materials for the manuscript:

**Longitudinal UAV remote sensing for wheat breeding decisions: cross-experiment transferability and early selection utility**

Current target journal: *Remote Sensing Applications: Society and Environment* (RSASE).

> Repository-path note: this directory keeps the original `CEA_UAV_phenomics_2026` path so previously shared links remain stable. The contents are now harmonized for the RSASE revision.

## Contents

- `data/CEA_wheat_phenomics_154_genotypes.csv` — genotype-level analytical dataset: 154 wheat genotypes, DH, GY, and 120 longitudinal spectral predictors from 10 UAV flights.
- `scripts/RSASE_master_pipeline.py` — production reference workflow for the audited RSASE revision.
- `scripts/CEA_master_pipeline.R` — companion R implementation retained for transparency; it is not the production source of the reported numerical results.

## Experimental domains

- E1: 49 cultivars
- E2: 56 breeding populations
- E3: 19 breeding populations
- E4: 30 elite lines

The experiments are treated as contrasting predictive domains. Differences among them should not be interpreted as isolated causal breeding-stage effects.

## Validation principles

1. Data-dependent preprocessing is fitted only within training resamples.
2. Within-experiment and cross-experiment prediction are distinct validation targets.
3. Cross-experiment transfer reports both **Pearson r** (linear concordance) and **Spearman rho** (ranking transfer), together with R2, RMSE, and MAE.
4. Longitudinal prediction uses cumulative UAV windows F1-F10.
5. The longitudinal XGBoost specification is frozen for the RSASE revision to isolate temporal information accumulation from repeated hyperparameter re-optimization.

## Frozen longitudinal XGBoost specification

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

A targeted audit with this fixed specification preserved the substantive longitudinal conclusions: E1-GY remains actionable at F6, E1-DH remains late, and weak E3/E4 signals remain weak.

## Related prior work

The E1 field campaign overlaps with the 49-genotype longitudinal dataset previously analyzed in:

Silva, C. M. et al. (2026). *Optimizing indirect selection of tropical wheat genotypes using high-throughput longitudinal phenotyping and trait relationships*. **Remote Sensing Applications: Society and Environment**, 41, 101892. https://doi.org/10.1016/j.rsase.2026.101892

That paper focused on longitudinal repeatability, trait relationships, and indirect selection. The present study addresses distinct questions of leakage-safe prediction, cross-experiment transferability, temporal information accumulation, and selection-decision stability.

## Data and code availability

The processed analytical dataset and the audited analysis workflow are publicly available in this repository. Additional materials can be obtained from the corresponding author upon reasonable request.

Corresponding author: **Maicon Nardino**, Federal University of Viçosa, Brazil.
