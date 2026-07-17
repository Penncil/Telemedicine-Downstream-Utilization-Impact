library(data.table)
library(dplyr)
library(xgboost)
library(ggplot2)
library(tableone)
library(openxlsx)

horizon <- "3m"
stopifnot(horizon %in% c("3m", "6m"))
followup_label <- if (horizon == "3m") "3-month" else "6-month"

data_dir    <- "PATH/TO/DATA"
results_dir <- "PATH/TO/RESULTS"

tele_list <- c("Telehealth", "Telemedicine", "Telemedicine (Non-Chargeable)",
               "Virtual Visit")
high_virtual_domains <- c("END", "MBD", "SKN", "MUS")

covars <- c("AGE_ENC", "GENDER", "race_ethnicity", "fin_class_category",
            "CCI", "n_inpatient_prior_cat", "n_ed_prior_cat", "distance",
            "CCSR_broad", "cohort_entry_month")

outcome_total   <- "n_followup_same_ccsr"
outcome_tele_fu <- "n_followup_tele"
outcome_f2f_fu  <- "n_followup_f2f"

xgb_params <- list(objective = "count:poisson", eval_metric = "poisson-nloglik",
                   max_depth = 4, eta = 0.1, subsample = 0.8,
                   colsample_bytree = 0.8, min_child_weight = 10)
nrounds <- 500
early_stopping <- 20
seed  <- 42
delta <- 0.25

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
wb <- createWorkbook()

in_path <- file.path(data_dir, paste0("telemedicine_cohort_enriched_", horizon, ".csv"))
stopifnot(file.exists(in_path))
df <- fread(in_path)

df[, treatment  := fifelse(ENCOUNTER_TYPE %in% tele_list, 1L, 0L)]
df[, CCSR_broad := substr(INDEX_CCSR, 1, 3)]

n_before_virtual <- nrow(df)
df <- df[CCSR_broad %in% high_virtual_domains]
cat("Virtualizability filter: kept", nrow(df), "of", n_before_virtual,
    "(", round(nrow(df) / n_before_virtual * 100, 1), "%)\n")

df[, age_group := cut(AGE_ENC, breaks = c(18, 35, 50, 65, Inf),
                      labels = c("18-34", "35-49", "50-64", "65+"), right = FALSE)]
df[, n_ed_prior_cat := fcase(n_ed_prior == 0, "0", n_ed_prior == 1, "1",
                             n_ed_prior == 2, "2", n_ed_prior >= 3, "2+")]
df[, n_inpatient_prior_cat := fcase(n_inpatient_prior == 0, "0", n_inpatient_prior == 1, "1",
                                    n_inpatient_prior == 2, "2", n_inpatient_prior >= 3, "2+")]

n_before_complete <- nrow(df)
df <- df[complete.cases(df[, ..covars])]
cat("Complete-case filter: kept", nrow(df), "of", n_before_complete,
    "(dropped", n_before_complete - nrow(df), ")\n")

table1_vars <- c("AGE_ENC", "age_group", "GENDER", "race_ethnicity", "fin_class_category",
                 "CCI", "n_inpatient_prior_cat", "n_ed_prior_cat", "distance", "CCSR_broad",
                 "n_followup_same_ccsr", "n_followup_tele", "n_followup_f2f")
cat_vars <- c("age_group", "GENDER", "race_ethnicity", "fin_class_category",
              "CCSR_broad", "n_inpatient_prior_cat", "n_ed_prior_cat")

df[, treatment_label := fifelse(treatment == 1L, "Telemedicine", "F2F")]
tab1 <- CreateTableOne(vars = table1_vars, strata = "treatment_label",
                       data = as.data.frame(df), factorVars = cat_vars, test = TRUE)
tab1_print <- print(tab1, showAllLevels = TRUE, smd = TRUE, printToggle = FALSE)
addWorksheet(wb, "Table1")
writeData(wb, "Table1", as.data.frame(tab1_print), rowNames = TRUE)

prep_xgb_data <- function(data, covars, outcome) {
  dt <- as.data.frame(data)[, c(covars, outcome), drop = FALSE]
  complete <- complete.cases(dt)
  dt <- dt[complete, , drop = FALSE]
  mm <- model.matrix(as.formula(paste0("~ ", paste(covars, collapse = " + "), " - 1")),
                     data = dt)
  list(dmat = xgb.DMatrix(data = mm, label = dt[[outcome]]),
       matrix = mm, complete_idx = which(complete))
}

prep_all    <- prep_xgb_data(df, covars, outcome_total)
df_complete <- df[prep_all$complete_idx]
mm_all      <- prep_all$matrix
label_all   <- df_complete[[outcome_total]]
dmat_all    <- prep_all$dmat

set.seed(seed)
tele_rows      <- which(df_complete$treatment == 1L)
f2f_rows       <- which(df_complete$treatment == 0L)
tele_train_idx <- sample(tele_rows, 0.8 * length(tele_rows))
f2f_train_idx  <- sample(f2f_rows,  0.8 * length(f2f_rows))
tele_test_idx  <- setdiff(tele_rows, tele_train_idx)
f2f_test_idx   <- setdiff(f2f_rows,  f2f_train_idx)

make_dmat <- function(idx, outcome) {
  xgb.DMatrix(data = mm_all[idx, , drop = FALSE], label = df_complete[[outcome]][idx])
}

train_arm <- function(train_idx, test_idx, outcome, verbose = 0) {
  xgb.train(params = xgb_params, nrounds = nrounds, verbose = verbose,
            data = make_dmat(train_idx, outcome),
            watchlist = list(train = make_dmat(train_idx, outcome),
                             test  = make_dmat(test_idx,  outcome)),
            early_stopping_rounds = early_stopping)
}

xgb_tele <- train_arm(tele_train_idx, tele_test_idx, outcome_total, verbose = 1)
xgb_f2f  <- train_arm(f2f_train_idx,  f2f_test_idx,  outcome_total, verbose = 1)

xgb_tele_fu_tele <- train_arm(tele_train_idx, tele_test_idx, outcome_tele_fu)
xgb_tele_fu_f2f  <- train_arm(f2f_train_idx,  f2f_test_idx,  outcome_tele_fu)
xgb_f2f_fu_tele  <- train_arm(tele_train_idx, tele_test_idx, outcome_f2f_fu)
xgb_f2f_fu_f2f   <- train_arm(f2f_train_idx,  f2f_test_idx,  outcome_f2f_fu)

tele_train_dt <- df_complete[tele_train_idx]
tele_test_dt  <- df_complete[tele_test_idx]
f2f_train_dt  <- df_complete[f2f_train_idx]
f2f_test_dt   <- df_complete[f2f_test_idx]

glm_formula <- as.formula(paste0(outcome_total, " ~ ", paste(covars, collapse = " + ")))
glm_tele <- glm(glm_formula, data = tele_train_dt, family = poisson())
glm_f2f  <- glm(glm_formula, data = f2f_train_dt,  family = poisson())

eval_model <- function(observed, predicted, label) {
  data.table(model = label,
             MAE  = mean(abs(observed - predicted)),
             RMSE = sqrt(mean((observed - predicted)^2)),
             mean_observed  = mean(observed),
             mean_predicted = mean(predicted))
}

pred_xgb_tele_test <- predict(xgb_tele, make_dmat(tele_test_idx, outcome_total))
pred_xgb_f2f_test  <- predict(xgb_f2f,  make_dmat(f2f_test_idx,  outcome_total))
pred_glm_tele_test <- predict(glm_tele, newdata = tele_test_dt, type = "response")
pred_glm_f2f_test  <- predict(glm_f2f,  newdata = f2f_test_dt,  type = "response")
obs_tele_test <- label_all[tele_test_idx]
obs_f2f_test  <- label_all[f2f_test_idx]

eval_results <- rbindlist(list(
  eval_model(obs_tele_test, pred_xgb_tele_test, "XGBoost-Tele"),
  eval_model(obs_f2f_test,  pred_xgb_f2f_test,  "XGBoost-F2F"),
  eval_model(obs_tele_test, pred_glm_tele_test, "GLM-Tele"),
  eval_model(obs_f2f_test,  pred_glm_f2f_test,  "GLM-F2F")
))
print(eval_results)
addWorksheet(wb, "ModelEval")
writeData(wb, "ModelEval", as.data.frame(eval_results))

make_calibration <- function(observed, predicted, group_label) {
  dt <- data.table(obs = observed, pred = predicted)
  dt[, decile := ntile(pred, 10)]
  dt[, .(mean_obs = mean(obs), mean_pred = mean(pred), group = group_label), by = decile]
}
cal_data <- rbindlist(list(
  make_calibration(obs_tele_test, pred_xgb_tele_test, "Tele"),
  make_calibration(obs_f2f_test,  pred_xgb_f2f_test,  "F2F")
))
p_cal <- ggplot(cal_data, aes(x = mean_pred, y = mean_obs, color = group)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = paste0("XGBoost Calibration Plot (", followup_label, ", Test Set)"),
       x = "Mean Predicted", y = "Mean Observed", color = "Group") +
  theme_minimal()
ggsave(file.path(results_dir, paste0("calibration_plot_", horizon, ".png")),
       p_cal, width = 7, height = 5)

df_complete[, pred_tele := predict(xgb_tele, dmat_all)]
df_complete[, pred_f2f  := predict(xgb_f2f,  dmat_all)]
df_complete[, elastic_util := pred_tele - pred_f2f]

df_complete[, pred_tele_fu_tele := predict(xgb_tele_fu_tele, dmat_all)]
df_complete[, pred_tele_fu_f2f  := predict(xgb_tele_fu_f2f,  dmat_all)]
df_complete[, pred_f2f_fu_tele  := predict(xgb_f2f_fu_tele,  dmat_all)]
df_complete[, pred_f2f_fu_f2f   := predict(xgb_f2f_fu_f2f,   dmat_all)]

df_complete[, eu_tele_channel := pred_tele_fu_tele - pred_tele_fu_f2f]
df_complete[, eu_f2f_channel  := pred_f2f_fu_tele  - pred_f2f_fu_f2f]

overall_summary <- data.table(
  metric = c("N total", "N complete", "N tele", "N f2f",
             "Mean EU (total)", "Median EU (total)", "SD EU (total)",
             "Mean EU (tele channel)", "Mean EU (f2f channel)",
             "Mean EU (tele + f2f)", "Additivity gap",
             "Pct EU > 0", "Pct EU > 0.3"),
  value = c(nrow(df), nrow(df_complete),
            sum(df_complete$treatment == 1L), sum(df_complete$treatment == 0L),
            round(mean(df_complete$elastic_util), 4),
            round(median(df_complete$elastic_util), 4),
            round(sd(df_complete$elastic_util), 4),
            round(mean(df_complete$eu_tele_channel), 4),
            round(mean(df_complete$eu_f2f_channel), 4),
            round(mean(df_complete$eu_tele_channel) + mean(df_complete$eu_f2f_channel), 4),
            round(mean(df_complete$elastic_util) -
                    (mean(df_complete$eu_tele_channel) + mean(df_complete$eu_f2f_channel)), 4),
            round(mean(df_complete$elastic_util > 0) * 100, 1),
            round(mean(df_complete$elastic_util > 0.3) * 100, 1))
)
print(overall_summary)
addWorksheet(wb, "OverallSummary")
writeData(wb, "OverallSummary", as.data.frame(overall_summary))

p_elastic <- ggplot(df_complete, aes(x = elastic_util)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = paste0("Distribution of Downstream Utilization Impact (",
                      followup_label, ", High Virtualizability)"),
       x = "Predicted DUI (Telemedicine - F2F)", y = "Count") +
  theme_minimal()
ggsave(file.path(results_dir, paste0("elastic_util_dist_", horizon, ".png")),
       p_elastic, width = 8, height = 5)

hetero_ccsr <- df_complete[, .(
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
print(hetero_ccsr)
addWorksheet(wb, "Heterogeneity_CCSR")
writeData(wb, "Heterogeneity_CCSR", as.data.frame(hetero_ccsr))

p_hetero <- ggplot(hetero_ccsr, aes(x = reorder(CCSR_broad, mean_eu), y = mean_eu)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Downstream Utilization Impact by CCSR Domain (High Virtualizability)",
       x = "CCSR Domain", y = "Mean DUI") +
  theme_minimal()
ggsave(file.path(results_dir, paste0("heterogeneity_ccsr_plot_", horizon, ".png")),
       p_hetero, width = 8, height = 5)

decomp_long <- melt(
  hetero_ccsr[, .(CCSR_broad, `Tele Channel` = mean_eu_tele, `F2F Channel` = mean_eu_f2f)],
  id.vars = "CCSR_broad", variable.name = "Channel", value.name = "Effect"
)
p_channel <- ggplot(decomp_long, aes(x = reorder(CCSR_broad, -Effect), y = Effect, fill = Channel)) +
  geom_col(position = "stack") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Tele Channel" = "steelblue", "F2F Channel" = "coral")) +
  labs(title = "DUI Decomposition by Channel and CCSR Domain",
       x = "CCSR Domain", y = paste0("Mean Effect (visits per ", followup_label, ")")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(file.path(results_dir, paste0("heterogeneity_channel_decomp_", horizon, ".png")),
       p_channel, width = 8, height = 5)

df_tele_2x2 <- df_complete[treatment == 1L]
n_before_washout <- nrow(df_tele_2x2)
df_tele_2x2 <- df_tele_2x2[!(as.logical(prior_same_ccsr_has_tele) == TRUE)]
cat("Washout removed:", n_before_washout - nrow(df_tele_2x2), "\n")

df_2x2 <- df_tele_2x2[!is.na(first_followup_modality)]
df_2x2[, behavioral_quadrant := fcase(
  n_followup_tele >= 1 & n_followup_f2f >= 1, "Q2: Both",
  n_followup_tele >= 1 & n_followup_f2f == 0, "Q1: Tele only",
  n_followup_tele == 0 & n_followup_f2f >= 1, "Q3: F2F only",
  n_followup_tele == 0 & n_followup_f2f == 0, "Q4: None"
)]
cat("Patients in behavioral 2x2:", nrow(df_2x2), "\n")

quad_summary <- df_2x2[, .(n = .N, pct = round(.N / nrow(df_2x2) * 100, 1)),
                       by = behavioral_quadrant][order(behavioral_quadrant)]
print(quad_summary)
addWorksheet(wb, "Behavioral2x2_Summary")
writeData(wb, "Behavioral2x2_Summary", as.data.frame(quad_summary))

quad_ccsr <- df_2x2[, .(n = .N), by = .(behavioral_quadrant, CCSR_broad)]
quad_ccsr[, pct := round(n / sum(n) * 100, 1), by = behavioral_quadrant]
quad_ccsr <- quad_ccsr[order(behavioral_quadrant, -n)]
addWorksheet(wb, "Behavioral2x2_CCSR")
writeData(wb, "Behavioral2x2_CCSR", as.data.frame(quad_ccsr))

df_complete[, decomp_quadrant := fcase(
  eu_tele_channel >  delta & eu_f2f_channel <= delta, "Q1: Convenience-driven",
  eu_tele_channel >  delta & eu_f2f_channel >  delta, "Q2: Both channels up",
  eu_tele_channel <= delta & eu_f2f_channel >  delta, "Q3: Escalation-dominant",
  eu_tele_channel <= delta & eu_f2f_channel <= delta, "Q4: No meaningful effect"
)]

decomp_summary <- df_complete[, .(
  n = .N, pct = round(.N / nrow(df_complete) * 100, 1),
  mean_eu_tele = round(mean(eu_tele_channel), 4),
  mean_eu_f2f  = round(mean(eu_f2f_channel), 4),
  mean_eu      = round(mean(elastic_util), 4)
), by = decomp_quadrant][order(decomp_quadrant)]
print(decomp_summary)

decomp_ccsr <- df_complete[, .(n = .N), by = .(decomp_quadrant, CCSR_broad)]
decomp_ccsr[, pct := round(n / sum(n) * 100, 1), by = decomp_quadrant]
decomp_ccsr <- decomp_ccsr[order(decomp_quadrant, -n)]

addWorksheet(wb, "Decomp2x2_Summary")
writeData(wb, "Decomp2x2_Summary", as.data.frame(decomp_summary))
addWorksheet(wb, "Decomp2x2_CCSR")
writeData(wb, "Decomp2x2_CCSR", as.data.frame(decomp_ccsr))

xlsx_path <- file.path(results_dir, paste0("main_results_", horizon, ".xlsx"))
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat("Saved:", xlsx_path, "\n")

fwrite(df_complete, file.path(data_dir, paste0("telemedicine_analysis_final_", horizon, ".csv")))
cat("Saved: telemedicine_analysis_final_", horizon, ".csv\n", sep = "")
