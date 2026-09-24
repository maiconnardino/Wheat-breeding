# ==============================================================================
# CEA UAV PHENOMICS 2026 — REPRODUCIBILITY WORKFLOW
# Manuscript: Longitudinal UAV phenomics for wheat breeding decisions:
#             cross-experiment transferability and early selection utility
# Target journal: Computers and Electronics in Agriculture
#
# This script is designed to accompany:
# data/CEA_wheat_phenomics_154_genotypes.csv
#
# PRINCIPLES
# 1. All data-dependent preprocessing is estimated from training data only.
# 2. Within-experiment and cross-experiment validation are separate targets.
# 3. Longitudinal prediction uses cumulative flight windows F1...F10.
# 4. Breeding utility is evaluated by ranking and Top-k recovery.
# ==============================================================================

# 00. Packages ----------------------------------------------------------------
pkgs <- c(
  "readr","dplyr","tidyr","purrr","stringr","tibble",
  "rsample","recipes","parsnip","workflows","tune","dials","yardstick",
  "glmnet","ranger","kernlab","xgboost","ggplot2"
)

missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(missing_pkgs) > 0) install.packages(missing_pkgs, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(20260922)

# 01. Data --------------------------------------------------------------------
DATA_FILE <- file.path("data", "CEA_wheat_phenomics_154_genotypes.csv")

dat <- readr::read_csv(DATA_FILE, show_col_types = FALSE)

stopifnot(nrow(dat) == 154)
stopifnot(dplyr::n_distinct(dat$GEN) == 154)
stopifnot(all(c("GEN","Genetic material","DH","GY") %in% names(dat)))

# Reconstruct experiment membership from the archived anonymized GEN mapping:
# E1 = G1-G49; E2 = G50-G105; E3 = G106-G124; E4 = G125-G154.
dat <- dat |>
  dplyr::mutate(
    GEN_num = as.integer(stringr::str_remove(GEN, "^G")),
    Stage = dplyr::case_when(
      GEN_num <= 49  ~ "E1",
      GEN_num <= 105 ~ "E2",
      GEN_num <= 124 ~ "E3",
      TRUE           ~ "E4"
    ),
    Stage = factor(Stage, levels = c("E1","E2","E3","E4"))
  )

stopifnot(all(table(dat$Stage) == c(E1=49,E2=56,E3=19,E4=30)))

phenomic_cols <- names(dat)[stringr::str_detect(names(dat), "flight[0-9]{2}$")]
stopifnot(length(phenomic_cols) == 120)

flight_num <- function(x){
  as.integer(stringr::str_extract(x, "(?<=flight)[0-9]{2}$"))
}

# 02. Metrics -----------------------------------------------------------------
pearson_r <- function(truth, estimate){
  suppressWarnings(stats::cor(truth, estimate, use = "complete.obs"))
}

metrics_from_predictions <- function(truth, pred){
  tibble::tibble(
    n = length(truth),
    r = pearson_r(truth, pred),
    R2 = 1 - sum((truth-pred)^2) / sum((truth-mean(truth))^2),
    RMSE = sqrt(mean((truth-pred)^2)),
    MAE = mean(abs(truth-pred))
  )
}

# 03. Leakage-safe recipe ------------------------------------------------------
make_recipe <- function(train_data, outcome, predictors){
  f <- stats::reformulate(predictors, response = outcome)
  recipes::recipe(f, data = train_data) |>
    recipes::step_zv(recipes::all_predictors()) |>
    recipes::step_corr(recipes::all_predictors(), threshold = 0.99) |>
    recipes::step_normalize(recipes::all_predictors())
}

# 04. Model specifications -----------------------------------------------------
ridge_spec <- parsnip::linear_reg(penalty = tune::tune(), mixture = 0) |>
  parsnip::set_engine("glmnet")

enet_spec <- parsnip::linear_reg(
  penalty = tune::tune(),
  mixture = tune::tune()
) |>
  parsnip::set_engine("glmnet")

rf_spec <- parsnip::rand_forest(
  mtry = tune::tune(),
  min_n = tune::tune(),
  trees = 1000
) |>
  parsnip::set_engine("ranger") |>
  parsnip::set_mode("regression")

svm_spec <- parsnip::svm_rbf(
  cost = tune::tune(),
  rbf_sigma = tune::tune(),
  margin = tune::tune()
) |>
  parsnip::set_engine("kernlab") |>
  parsnip::set_mode("regression")

xgb_spec <- parsnip::boost_tree(
  trees = tune::tune(),
  tree_depth = tune::tune(),
  learn_rate = tune::tune(),
  loss_reduction = tune::tune(),
  min_n = tune::tune(),
  sample_size = tune::tune(),
  mtry = tune::tune()
) |>
  parsnip::set_engine("xgboost", nthread = 1) |>
  parsnip::set_mode("regression")

model_specs <- list(
  Ridge = ridge_spec,
  ElasticNet = enet_spec,
  RF = rf_spec,
  SVM = svm_spec,
  XGBoost = xgb_spec
)

# 05. Within-experiment nested CV ---------------------------------------------
within_experiment_cv <- function(df, stage, outcome, flight_max = 10,
                                 outer_repeats = 3, seed = 20260922){

  set.seed(seed)
  d <- df |> dplyr::filter(Stage == stage)

  preds <- phenomic_cols[flight_num(phenomic_cols) <= flight_max]
  v_outer <- ifelse(nrow(d) < 40, 4, 5)

  outer <- rsample::vfold_cv(d, v = v_outer, repeats = outer_repeats)

  purrr::imap_dfr(outer$splits, function(split, fold_id){
    tr <- rsample::analysis(split)
    te <- rsample::assessment(split)

    rec <- make_recipe(tr, outcome, preds)
    inner <- rsample::vfold_cv(tr, v = min(4, max(2, floor(nrow(tr)/8))))

    purrr::imap_dfr(model_specs, function(spec, model_name){
      wf <- workflows::workflow() |>
        workflows::add_recipe(rec) |>
        workflows::add_model(spec)

      grid_size <- if(model_name %in% c("RF","XGBoost")) 20 else 15

      tuned <- tune::tune_grid(
        wf,
        resamples = inner,
        grid = grid_size,
        metrics = yardstick::metric_set(yardstick::rmse),
        control = tune::control_grid(verbose = FALSE)
      )

      best <- tune::select_best(tuned, metric = "rmse")
      final_wf <- tune::finalize_workflow(wf, best)
      fit_obj <- parsnip::fit(final_wf, data = tr)

      p <- predict(fit_obj, new_data = te)$.pred
      m <- metrics_from_predictions(te[[outcome]], p)

      dplyr::bind_cols(
        tibble::tibble(
          Stage = stage,
          Trait = outcome,
          Model = model_name,
          flight_max = flight_max,
          fold = fold_id
        ),
        m
      )
    })
  })
}

# 06. Cross-experiment transfer -----------------------------------------------
cross_experiment_transfer <- function(df, target_stage, outcome,
                                      flight_max = 10, seed = 20260922){

  set.seed(seed)

  train_df <- df |> dplyr::filter(Stage != target_stage)
  test_df  <- df |> dplyr::filter(Stage == target_stage)

  preds <- phenomic_cols[flight_num(phenomic_cols) <= flight_max]
  rec <- make_recipe(train_df, outcome, preds)

  # Tuning uses grouped folds defined by source experiment.
  inner <- rsample::group_vfold_cv(train_df, group = Stage, v = 3)

  purrr::imap_dfr(model_specs, function(spec, model_name){
    wf <- workflows::workflow() |>
      workflows::add_recipe(rec) |>
      workflows::add_model(spec)

    grid_size <- if(model_name %in% c("RF","XGBoost")) 25 else 20

    tuned <- tune::tune_grid(
      wf,
      resamples = inner,
      grid = grid_size,
      metrics = yardstick::metric_set(yardstick::rmse),
      control = tune::control_grid(verbose = FALSE)
    )

    best <- tune::select_best(tuned, metric = "rmse")
    final_wf <- tune::finalize_workflow(wf, best)
    fit_obj <- parsnip::fit(final_wf, data = train_df)

    p <- predict(fit_obj, new_data = test_df)$.pred
    m <- metrics_from_predictions(test_df[[outcome]], p)

    dplyr::bind_cols(
      tibble::tibble(
        TargetStage = target_stage,
        Trait = outcome,
        Model = model_name,
        flight_max = flight_max
      ),
      m
    )
  })
}

# 07. Longitudinal cumulative windows -----------------------------------------
run_longitudinal <- function(df, outer_repeats = 3){
  grid <- tidyr::crossing(
    Stage = levels(df$Stage),
    Trait = c("DH","GY"),
    flight_max = 1:10
  )

  purrr::pmap_dfr(
    grid,
    function(Stage,Trait,flight_max){
      within_experiment_cv(
        df, stage = Stage, outcome = Trait,
        flight_max = flight_max,
        outer_repeats = outer_repeats,
        seed = 20260922 + flight_max
      )
    }
  )
}

# 08. Selection utility --------------------------------------------------------
topk_recovery <- function(truth, pred, p = 0.20){
  n <- length(truth)
  k <- max(1, ceiling(n*p))
  obs <- order(truth, decreasing = TRUE)[seq_len(k)]
  est <- order(pred, decreasing = TRUE)[seq_len(k)]

  tibble::tibble(
    selection_fraction = p,
    k = k,
    recovery = length(intersect(obs,est))/k,
    random_expectation = k/n,
    enrichment = (length(intersect(obs,est))/k)/(k/n)
  )
}

# 09. Example execution --------------------------------------------------------
# These calls may be computationally intensive.
#
# within_all <- tidyr::crossing(
#   Stage = c("E1","E2","E3","E4"),
#   Trait = c("DH","GY")
# ) |>
#   purrr::pmap_dfr(~within_experiment_cv(dat, ..1, ..2))
#
# cross_all <- tidyr::crossing(
#   TargetStage = c("E1","E2","E3","E4"),
#   Trait = c("DH","GY")
# ) |>
#   purrr::pmap_dfr(~cross_experiment_transfer(dat, ..1, ..2))
#
# longitudinal <- run_longitudinal(dat, outer_repeats = 3)
#
# readr::write_csv(within_all, "within_experiment_results.csv")
# readr::write_csv(cross_all, "cross_experiment_results.csv")
# readr::write_csv(longitudinal, "longitudinal_results.csv")

message("CEA UAV phenomics reproducibility workflow loaded successfully.")