# ==============================================================================
# COMPANION R WORKFLOW — RSASE REVISION
# ==============================================================================
# Manuscript:
# Longitudinal UAV remote sensing for wheat breeding decisions:
# cross-experiment transferability and early selection utility
#
# IMPORTANT:
# The production reference for the reported RSASE revision is:
#   scripts/RSASE_master_pipeline.py
#
# This R file is retained as a companion implementation for transparency and
# teaching/reuse. It must not be cited as the script that generated the exact
# numerical tables unless its output has first been reconciled with the Python
# production workflow.
#
# Fixed longitudinal XGBoost specification used in the RSASE revision:
# nrounds/trees = 120
# max_depth = 2
# eta/learn_rate = 0.05
# subsample = 0.80
# colsample_bytree = 0.80
# min_child_weight = 1
# lambda = 1.0
# alpha = 0.0
# gamma = 0.0
# seed = 20260922
#
# Cross-experiment reporting should include:
# Pearson r, Spearman rho, R2, RMSE and MAE.
# ==============================================================================

message(
  "This is the companion R workflow. Use scripts/RSASE_master_pipeline.py ",
  "as the production reference for the RSASE revision."
)
