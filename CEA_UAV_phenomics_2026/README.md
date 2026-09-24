# Longitudinal UAV phenomics for wheat breeding decisions

Reproducibility materials for the manuscript:

**Longitudinal UAV phenomics for wheat breeding decisions: cross-experiment transferability and early selection utility**

Target journal: *Computers and Electronics in Agriculture*.

## Contents

- `data/CEA_wheat_phenomics_154_genotypes.csv` — genotype-level analytical dataset used for phenomic prediction (154 wheat genotypes; DH, GY, and 120 longitudinal spectral predictors from 10 UAV flights).
- `scripts/CEA_master_pipeline.R` — R reproducibility workflow implementing leakage-safe preprocessing, within-experiment prediction, cross-experiment transfer, longitudinal prediction, and breeding-selection utility.

## Experimental domains

The analytical dataset comprises four contrasting experiments:

- E1: 49 cultivars
- E2: 56 breeding populations
- E3: 19 breeding populations
- E4: 30 elite lines

The experiment labels can be reconstructed from the anonymized genotype IDs as documented in the script.

## Reproducibility note

All data-dependent preprocessing and model tuning are designed to occur only within training resamples. The experiments are treated as contrasting predictive domains; differences among them should not be interpreted as isolated causal breeding-stage effects.

## Data and code availability

The processed data and analysis workflow are publicly available in this repository. Additional material can be obtained from the corresponding author upon reasonable request.

Corresponding author: **Maicon Nardino**, Federal University of Viçosa, Brazil.

## Citation

Please cite the associated article when using these materials.