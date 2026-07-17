library(data.table)
library(xgboost)
library(openxlsx)

data_dir    <- "PATH/TO/DATA"
results_dir <- "PATH/TO/RESULTS"
xwalk_path  <- "PATH/TO/Specialty_Categorization.csv"

horizons  <- c("3m", "6m")
MIN_ARM_N <- 200
DELTA_3M  <- 0.25
SEED      <- 42

covars <- c("AGE_ENC", "GENDER", "race_ethnicity", "fin_class_category",
            "CCI", "n_inpatient_prior_cat", "n_ed_prior_cat", "distance",
            "CCSR_broad", "cohort_entry_month")
outcome_total   <- "n_followup_same_ccsr"
outcome_tele_fu <- "n_followup_tele"
outcome_f2f_fu  <- "n_followup_f2f"

xgb_params <- list(objective = "count:poisson", eval_metric = "poisson-nloglik",
                   max_depth = 4, eta = 0.1, subsample = 0.8,
                   colsample_bytree = 0.8, min_child_weight = 10)

add_care_level <- function(dt, xwalk_path) {
  stopifnot(file.exists(xwalk_path))
  xw <- fread(xwalk_path)
  stopifnot(all(c("specialty", "care_category") %in% names(xw)))
  xw <- unique(xw[, .(specialty = trimws(specialty),
                      care_category = trimws(care_category))])
  dt[, spec_key := trimws(as.character(PROV_SPECIALTY))]
  dt <- merge(dt, xw, by.x = "spec_key", by.y = "specialty", all.x = TRUE)
  dt[, is_primary_care := fcase(care_category == "Primary Care",   1L,
                                care_category == "Specialty Care", 0L,
                                default = NA_integer_)]
  dt[, spec_key := NULL]
  dt[]
}

unmatched_report <- function(dt) {
  dt[is.na(care_category), .(n = .N), by = .(PROV_SPECIALTY)][order(-n)]
}

get_complete <- function(df_sub) {
  d  <- as.data.frame(df_sub)[, c(covars, outcome_total), drop = FALSE]
  cc <- complete.cases(d)
  d  <- d[cc, , drop = FALSE]
  mm <- model.matrix(
    as.formula(paste0("~ ", paste(covars, collapse = " + "), " - 1")), data = d)
  list(dfc = df_sub[which(cc)], mm = mm, y = d[[outcome_total]])
}

fit_scores <- function(dfc, mm, y, seed = SEED) {
  tele_rows <- which(dfc$treatment == 1L)
  f2f_rows  <- which(dfc$treatment == 0L)

  set.seed(seed)
  tele_tr <- sample(tele_rows, floor(0.8 * length(tele_rows)))
  f2f_tr  <- sample(f2f_rows,  floor(0.8 * length(f2f_rows)))
  tele_te <- setdiff(tele_rows, tele_tr)
  f2f_te  <- setdiff(f2f_rows,  f2f_tr)

  dm <- function(idx, oc) xgb.DMatrix(mm[idx, , drop = FALSE], label = dfc[[oc]][idx])
  tp <- function(tr, te, oc) {
    xgb.train(params = xgb_params, nrounds = 500, verbose = 0,
              data = dm(tr, oc),
              watchlist = list(train = dm(tr, oc), test = dm(te, oc)),
              early_stopping_rounds = 20)
  }

  m_tele <- tp(tele_tr, tele_te, outcome_total)
  m_f2f  <- tp(f2f_tr,  f2f_te,  outcome_total)
  m_tt   <- tp(tele_tr, tele_te, outcome_tele_fu)
  m_tf   <- tp(f2f_tr,  f2f_te,  outcome_tele_fu)
  m_ft   <- tp(tele_tr, tele_te, outcome_f2f_fu)
  m_ff   <- tp(f2f_tr,  f2f_te,  outcome_f2f_fu)

  dmat_all <- xgb.DMatrix(mm)
  dfc <- copy(dfc)
  dfc[, pred_tele       := predict(m_tele, dmat_all)]
  dfc[, pred_f2f        := predict(m_f2f,  dmat_all)]
  dfc[, elastic_util    := pred_tele - pred_f2f]
  dfc[, eu_tele_channel := predict(m_tt, dmat_all) - predict(m_tf, dmat_all)]
  dfc[, eu_f2f_channel  := predict(m_ft, dmat_all) - predict(m_ff, dmat_all)]

  ev <- function(obs, pred, l) {
    data.table(model     = l,
               MAE       = round(mean(abs(obs - pred)), 4),
               RMSE      = round(sqrt(mean((obs - pred)^2)), 4),
               mean_obs  = round(mean(obs), 4),
               mean_pred = round(mean(pred), 4))
  }
  eval_tab <- rbindlist(list(
    ev(y[tele_te], predict(m_tele, xgb.DMatrix(mm[tele_te, , drop = FALSE])), "XGB-Tele"),
    ev(y[f2f_te],  predict(m_f2f,  xgb.DMatrix(mm[f2f_te, , drop = FALSE])),  "XGB-F2F")))

  list(scored = dfc, eval = eval_tab)
}

overall_summary <- function(dfc) {
  data.table(
    metric = c("N complete", "N tele", "N f2f",
               "Mean EU (total)", "Median EU (total)", "SD EU (total)",
               "Mean EU (tele channel)", "Mean EU (f2f channel)",
               "Mean EU (tele + f2f)", "Additivity gap",
               "Pct EU > 0", "Pct EU > 0.3"),
    value = c(nrow(dfc), sum(dfc$treatment == 1L), sum(dfc$treatment == 0L),
              round(mean(dfc$elastic_util), 4),
              round(median(dfc$elastic_util), 4),
              round(sd(dfc$elastic_util), 4),
              round(mean(dfc$eu_tele_channel), 4),
              round(mean(dfc$eu_f2f_channel), 4),
              round(mean(dfc$eu_tele_channel) + mean(dfc$eu_f2f_channel), 4),
              round(mean(dfc$elastic_util) -
                      (mean(dfc$eu_tele_channel) + mean(dfc$eu_f2f_channel)), 4),
              round(mean(dfc$elastic_util > 0) * 100, 1),
              round(mean(dfc$elastic_util > 0.3) * 100, 1)))
}

heterogeneity_ccsr <- function(dfc) {
  dfc[, .(
    n              = .N,
    pct_tele       = round(mean(treatment == 1L) * 100, 1),
    mean_eu        = round(mean(elastic_util), 4),
    se_eu          = round(sd(elastic_util) / sqrt(.N), 4),
    ci_lower       = round(mean(elastic_util) - 1.96 * sd(elastic_util) / sqrt(.N), 4),
    ci_upper       = round(mean(elastic_util) + 1.96 * sd(elastic_util) / sqrt(.N), 4),
    mean_eu_tele   = round(mean(eu_tele_channel), 4),
    se_eu_tele     = round(sd(eu_tele_channel) / sqrt(.N), 4),
    mean_eu_f2f    = round(mean(eu_f2f_channel), 4),
    se_eu_f2f      = round(sd(eu_f2f_channel) / sqrt(.N), 4),
    sum_check      = round(mean(eu_tele_channel) + mean(eu_f2f_channel), 4),
    additivity_gap = round(mean(elastic_util) -
                             (mean(eu_tele_channel) + mean(eu_f2f_channel)), 4)
  ), by = CCSR_broad][order(-n)]
}

classify_decomp <- function(df, delta) {
  out <- copy(df)
  out[, decomp_quadrant := fcase(
    eu_tele_channel >  delta & eu_f2f_channel <= delta, "Q1: Convenience-driven",
    eu_tele_channel >  delta & eu_f2f_channel >  delta, "Q2: Both channels up",
    eu_tele_channel <= delta & eu_f2f_channel >  delta, "Q3: Escalation-dominant",
    eu_tele_channel <= delta & eu_f2f_channel <= delta, "Q4: No meaningful effect")]
  out
}

decomp_overall <- function(df, horizon, delta) {
  df[, .(horizon = horizon,
         delta   = round(delta, 4),
         n       = .N,
         pct     = round(.N / nrow(df) * 100, 2),
         mean_eu_tele = round(mean(eu_tele_channel), 4),
         mean_eu_f2f  = round(mean(eu_f2f_channel),  4),
         mean_eu      = round(mean(elastic_util),    4)),
     by = decomp_quadrant][order(decomp_quadrant)]
}

decomp_by_ccsr <- function(df, horizon, delta) {
  out <- df[, .(n = .N,
                mean_eu_tele = round(mean(eu_tele_channel), 4),
                mean_eu_f2f  = round(mean(eu_f2f_channel),  4),
                mean_eu      = round(mean(elastic_util),    4)),
            by = .(CCSR_broad, decomp_quadrant)]
  out[, total := sum(n), by = CCSR_broad]
  out[, pct   := round(n / total * 100, 2)]
  out[, total := NULL]
  out[, `:=`(horizon = horizon, delta = round(delta, 4))]
  setcolorder(out, c("horizon", "delta", "CCSR_broad", "decomp_quadrant",
                     "n", "pct", "mean_eu_tele", "mean_eu_f2f", "mean_eu"))
  out[order(CCSR_broad, decomp_quadrant)]
}

behavioral_2x2 <- function(dfc) {
  need <- c("prior_same_ccsr_has_tele", "first_followup_modality",
            "n_followup_tele", "n_followup_f2f")
  if (!all(need %in% names(dfc))) return(NULL)
  d <- dfc[treatment == 1L]
  d <- d[!(as.logical(prior_same_ccsr_has_tele) == TRUE)]
  d <- d[!is.na(first_followup_modality)]
  if (!nrow(d)) return(NULL)
  d[, behavioral_quadrant := fcase(
    n_followup_tele >= 1 & n_followup_f2f >= 1, "Q2: Both",
    n_followup_tele >= 1 & n_followup_f2f == 0, "Q1: Tele only",
    n_followup_tele == 0 & n_followup_f2f >= 1, "Q3: F2F only",
    n_followup_tele == 0 & n_followup_f2f == 0, "Q4: None")]
  s <- d[, .(n = .N, pct = round(.N / nrow(d) * 100, 1)),
         by = behavioral_quadrant][order(behavioral_quadrant)]
  c2 <- d[, .(n = .N), by = .(behavioral_quadrant, CCSR_broad)]
  c2[, pct := round(n / sum(n) * 100, 1), by = behavioral_quadrant]
  c2 <- c2[order(behavioral_quadrant, -n)]
  list(summary = s, ccsr = c2)
}

composition <- function(df_aug) {
  d <- df_aug[!is.na(care_category)]
  overall <- d[, .(n = .N,
                   n_tele = sum(treatment == 1L),
                   n_f2f  = sum(treatment == 0L)), by = care_category]
  overall[, pct := round(n / sum(n) * 100, 1)]
  overall[, pct_tele := round(n_tele / n * 100, 1)]
  bydom <- d[, .(n = .N,
                 n_tele = sum(treatment == 1L),
                 n_f2f  = sum(treatment == 0L)),
             by = .(CCSR_broad, care_category)]
  bydom[, pct_within_domain := round(n / sum(n) * 100, 1), by = CCSR_broad]
  bydom <- bydom[order(CCSR_broad, -n)]
  list(overall = overall, by_domain = bydom)
}

run_subgroup <- function(df_sub, lab, tag, horizon, delta) {
  cat("Subgroup:", lab, "| horizon:", horizon, "| delta:", round(delta, 4), "\n")
  comp   <- get_complete(df_sub)
  dfc    <- comp$dfc
  n_tele <- sum(dfc$treatment == 1L)
  n_f2f  <- sum(dfc$treatment == 0L)
  cat("  complete-case N:", nrow(dfc), "| tele:", n_tele, "| f2f:", n_f2f, "\n")

  behav <- behavioral_2x2(dfc)

  if (n_tele < MIN_ARM_N || n_f2f < MIN_ARM_N) {
    cat("  models skipped: an arm is below MIN_ARM_N =", MIN_ARM_N, "\n")
    return(list(tag = tag, lab = lab, ok = FALSE, n_tele = n_tele, n_f2f = n_f2f,
                eval = NULL, overall = NULL, hetero = NULL,
                decomp_sum = NULL, decomp_ccsr = NULL, behav = behav))
  }

  fit    <- fit_scores(dfc, comp$mm, comp$y)
  scored <- fit$scored
  ov     <- overall_summary(scored)
  het    <- heterogeneity_ccsr(scored)
  scored <- classify_decomp(scored, delta)
  dc_sum  <- decomp_overall(scored, horizon, delta)
  dc_ccsr <- decomp_by_ccsr(scored, horizon, delta)

  print(ov)
  print(dc_sum)

  list(tag = tag, lab = lab, ok = TRUE, n_tele = n_tele, n_f2f = n_f2f,
       eval = fit$eval, overall = ov, hetero = het,
       decomp_sum = dc_sum, decomp_ccsr = dc_ccsr, behav = behav)
}

add_sheet <- function(wb, name, x) {
  if (is.null(x) || !nrow(x)) return(invisible())
  addWorksheet(wb, name)
  writeData(wb, name, as.data.frame(x))
}

write_subgroup_sheets <- function(wb, res) {
  t <- res$tag
  add_sheet(wb, paste0(t, "_OverallSummary"), res$overall)
  add_sheet(wb, paste0(t, "_Hetero_CCSR"),    res$hetero)
  add_sheet(wb, paste0(t, "_ModelEval"),      res$eval)
  add_sheet(wb, paste0(t, "_Decomp2x2"),      res$decomp_sum)
  add_sheet(wb, paste0(t, "_Decomp2x2_CCSR"), res$decomp_ccsr)
  if (!is.null(res$behav)) {
    add_sheet(wb, paste0(t, "_Behav2x2"),      res$behav$summary)
    add_sheet(wb, paste0(t, "_Behav2x2_CCSR"), res$behav$ccsr)
  }
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

mean_total_3m <- mean(fread(file.path(data_dir, "telemedicine_analysis_final_3m.csv"))$elastic_util)
mean_total_6m <- mean(fread(file.path(data_dir, "telemedicine_analysis_final_6m.csv"))$elastic_util)
delta_6m_scaled  <- mean_total_6m / mean_total_3m * DELTA_3M
DELTA_BY_HORIZON <- c("3m" = DELTA_3M, "6m" = delta_6m_scaled)

cat("mean total DUI 3m =", round(mean_total_3m, 4),
    "| 6m =", round(mean_total_6m, 4), "\n")
cat("delta 3m =", round(DELTA_BY_HORIZON[["3m"]], 4),
    "| 6m =", round(DELTA_BY_HORIZON[["6m"]], 4), "\n")

for (h in horizons) {
  cat("HORIZON:", h, "\n")
  in_path <- file.path(data_dir, paste0("telemedicine_analysis_final_", h, ".csv"))
  stopifnot(file.exists(in_path))

  df <- fread(in_path)
  df[, treatment := as.integer(treatment)]
  if ("prior_same_ccsr_has_tele" %in% names(df)) {
    df[, prior_same_ccsr_has_tele := as.logical(prior_same_ccsr_has_tele)]
  }

  df <- add_care_level(df, xwalk_path)

  n_total <- nrow(df)
  n_pc    <- sum(df$is_primary_care == 1L, na.rm = TRUE)
  n_sp    <- sum(df$is_primary_care == 0L, na.rm = TRUE)
  n_na    <- sum(is.na(df$care_category))
  cat("N =", n_total, "| Primary Care =", n_pc, "| Specialist =", n_sp,
      "| unmatched =", n_na, "\n")

  unmatched <- unmatched_report(df)
  if (nrow(unmatched)) print(unmatched)

  wb   <- createWorkbook()
  comp <- composition(df)
  add_sheet(wb, "GroupComposition",        comp$overall)
  add_sheet(wb, "GroupComposition_Domain", comp$by_domain)
  add_sheet(wb, "Unmatched_Specialty",     unmatched)

  subgroups <- list(list(lab = "Primary Care", tag = "PC", flag = 1L),
                    list(lab = "Specialist",   tag = "SP", flag = 0L))
  log_rows <- list()
  for (sg in subgroups) {
    res <- run_subgroup(df[is_primary_care == sg$flag], sg$lab, sg$tag,
                        horizon = h, delta = DELTA_BY_HORIZON[[h]])
    write_subgroup_sheets(wb, res)
    log_rows[[sg$tag]] <- data.table(subgroup = res$lab, tag = res$tag,
                                     n_tele = res$n_tele, n_f2f = res$n_f2f,
                                     models_fit = res$ok)
  }
  add_sheet(wb, "AnalysisLog", rbindlist(log_rows))

  xlsx_path <- file.path(results_dir, paste0("subgroup_results_", h, ".xlsx"))
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  cat("Saved:", xlsx_path, "\n")
}
