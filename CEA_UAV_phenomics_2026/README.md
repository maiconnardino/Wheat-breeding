# Longitudinal UAV phenomics for wheat breeding decisions

Reproducibility materials for the manuscript:

**Longitudinal UAV phenomics for wheat breeding decisions: cross-experiment transferability and early selection utility**

Current target journal: *Remote Sensing Applications: Society and Environment*.

## Contents

- `data/CEA_wheat_phenomics_154_genotypes.csv` — genotype-level analytical dataset used for phenomic prediction (154 wheat genotypes; DH, GY, and 120 longitudinal spectral predictors from 10 UAV flights).
- `scripts/CEA_master_pipeline.R` — companion R implementation of the analysis framework. The reported numerical results were generated during the audited Python workflow; the repository is being harmonized so that the exact production scripts used for the submitted tables and figures are archived alongside the final submission.

## Experimental domains

The analytical dataset comprises four contrasting experiments:

- E1: 49 cultivars
- E2: 56 breeding populations
- E3: 19 breeding populations
- E4: 30 elite lines

The experiment labels can be reconstructed from the anonymized genotype IDs as documented in the script.

## Reproducibility note

All data-dependent preprocessing and model tuning are designed to occur only within training resamples. The experiments are treated as contrasting predictive domains; differences among them should not be interpreted as isolated causal breeding-stage effects. The E1 field campaign overlaps with the 49-genotype longitudinal dataset previously analyzed in Silva et al. (2026, Remote Sensing Applications: Society and Environment, 41, 101892; DOI: 10.1016/j.rsase.2026.101892) for repeatability, trait relationships, and indirect selection. The present study addresses distinct questions of leakage-safe prediction, cross-experiment transferability, temporal information accumulation, and selection-decision stability.

## Data and code availability

The processed data and analysis workflow are publicly available in this repository. Additional material can be obtained from the corresponding author upon reasonable request.

Corresponding author: **Maicon Nardino**, Federal University of Viçosa, Brazil.

## Citation

Please cite the associated article when using these materials.