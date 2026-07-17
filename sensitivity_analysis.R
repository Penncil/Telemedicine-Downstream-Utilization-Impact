library(data.table)
library(xgboost)
library(ranger)
library(lightgbm)
library(openxlsx)

data_dir    <- "PATH/TO/DATA"
results_dir <- "PATH/TO/RESULTS"

horizons <- c("3m", "6m")
seed     <- 42

covars <- c("AGE_ENC", "GENDER", "race_ethnicity", "fin_class_category",
            "CCI", "n_inpatient_prior_cat", "n_ed_prior_cat", "distance",
            "CCSR_broad", "cohort_entry_month")
outcome    <- "n_followup_same_ccsr"
factor_vars <- c("GENDER", "race_ethnicity", "fin_class_category",
                 "n_inpatient_prior_cat", "n_ed_prior_cat",
                 "CCSR_broad", "cohort_entry_month")
model_names <- c("xgboost", "glm_poisson", "ranger", "lightgbm")
cat_order   <- c("Total", "MBD", "END", "SKN", "MUS")

fit_xgboost <- function(mm_train, label_train, mm_full) {
  m <- xgb.train(
    params  = list(objective = "count:poisson", max_depth = 4, eta = 0.1,
                   subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 10),
    data    = xgb.DMatrix(mm_train, label = label_train),
    nrounds = 500, verbose = 0
  )
  predict(m, xgb.DMatrix(mm_full))
}

fit_glm_poisson <- function(df_train, df_full, formula_str) {
  m <- glm(as.formula(formula_str), data = df_train, family = poisson())
  predict(m, newdata = df_full, type = "response")
}

fit_ranger <- function(df_train, df_full, formula_str) {
  m <- ranger(formula       = as.formula(formula_str),
              data          = df_train,
              num.trees     = 500,
              min.node.size = 10,
              num.threads   = 4,
              verbose       = FALSE)
  predict(m, data = df_full)$predictions
}

fit_lightgbm <- function(mm_train, label_train, mm_full) {
  m <- lgb.train(
    params = list(objective = "poisson", learning_rate = 0.1, max_depth = 4,
                  feature_fraction = 0.8, bagging_fraction = 0.8, bagging_freq = 1,
                  min_data_in_leaf = 10, verbosity = -1),
    data    = lgb.Dataset(mm_train, label = label_train),
    nrounds = 500
  )
  predict(m, mm_full)
}

predict_arm <- function(model_name, mm_full, label, df_for_glm, formula_str, train_idx) {
  switch(model_name,
    xgboost     = fit_xgboost(mm_full[train_idx, ], label[train_idx], mm_full),
    glm_poisson = fit_glm_poisson(df_for_glm[train_idx], df_for_glm, formula_str),
    ranger      = fit_ranger(df_for_glm[train_idx], df_for_glm, formula_str),
    lightgbm    = fit_lightgbm(mm_full[train_idx, ], label[train_idx], mm_full)
  )
}

run_sensitivity <- function(horizon) {
  cat("Horizon:", horizon, "\n")
  in_path <- file.path(data_dir, paste0("telemedicine_analysis_final_", horizon, ".csv"))
  stopifnot(file.exists(in_path))
  df <- fread(in_path)

  formula_str <- paste0(outcome, " ~ ", paste(covars, collapse = " + "))
  mm_full <- model.matrix(
    as.formula(paste0("~ ", paste(covars, collapse = " + "), " - 1")),
    data = as.data.frame(df)[, covars, drop = FALSE]
  )
  label <- df[[outcome]]

  set.seed(seed)
  tele_rows      <- which(df$treatment == 1L)
  f2f_rows       <- which(df$treatment == 0L)
  tele_train_idx <- sample(tele_rows, 0.8 * length(tele_rows))
  f2f_train_idx  <- sample(f2f_rows,  0.8 * length(f2f_rows))

  df_for_glm <- copy(df)
  for (v in factor_vars) df_for_glm[[v]] <- as.factor(df_for_glm[[v]])

  results <- list()
  for (mn in model_names) {
    cat("  Model:", mn, "...\n")
    t0 <- Sys.time()
    pt <- predict_arm(mn, mm_full, label, df_for_glm, formula_str, tele_train_idx)
    pf <- predict_arm(mn, mm_full, label, df_for_glm, formula_str, f2f_train_idx)
    cat("    elapsed:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

    dui <- pt - pf

    overall <- data.table(horizon  = horizon,
                          model    = mn,
                          category = "Total",
                          n        = nrow(df),
                          mean_dui = mean(dui),
                          se_dui   = sd(dui) / sqrt(nrow(df)))

    df_tmp <- copy(df)
    df_tmp[, dui_tmp := dui]
    by_ccsr <- df_tmp[, .(horizon  = horizon,
                          model    = mn,
                          n        = .N,
                          mean_dui = mean(dui_tmp),
                          se_dui   = sd(dui_tmp) / sqrt(.N)),
                      by = .(category = CCSR_broad)]
    setcolorder(by_ccsr, c("horizon", "model", "category", "n", "mean_dui", "se_dui"))

    results[[mn]] <- rbind(overall, by_ccsr)
  }
  rbindlist(results)
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

all <- rbindlist(lapply(horizons, run_sensitivity))
all[, mean_dui := round(mean_dui, 4)]
all[, se_dui   := round(se_dui,   4)]

wide <- dcast(all, horizon + category ~ model, value.var = "mean_dui")
wide[, category := factor(category, levels = cat_order)]
wide <- wide[order(horizon, category)]

all[, category := factor(category, levels = cat_order)]
all <- all[order(horizon, model, category)]

print(wide)

wb <- createWorkbook()
addWorksheet(wb, "long"); writeData(wb, "long", as.data.frame(all))
addWorksheet(wb, "wide"); writeData(wb, "wide", as.data.frame(wide))
out_xlsx <- file.path(results_dir, "sensitivity_models.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat("Saved:", out_xlsx, "\n")
