# Downstream Utilization Impact (DUI) of Telemedicine

[![R](https://img.shields.io/badge/R-%E2%89%A5%204.2-276DC3?style=flat-square&logo=r&logoColor=white)](https://www.r-project.org/)
[![CCSR](https://img.shields.io/badge/AHRQ%20HCUP%20CCSR-v2023--1-8E6BB0?style=flat-square)](https://hcup-us.ahrq.gov/)
[![License: MIT](https://img.shields.io/badge/License-MIT-E8A33D?style=flat-square)](LICENSE)

Analysis code accompanying *Downstream Utilization Impact of Telemedicine in Outpatient Care* (Lu et al.).

This repository contains the statistical analysis code used to estimate the
Downstream Utilization Impact (DUI) of telemedicine versus face-to-face index
visits, decompose it into telemedicine-channel and face-to-face-channel
components, and run the sensitivity and subgroup analyses reported in the
manuscript and supplement.

## Contents

| File | Purpose |
|---|---|
| `analysis.R` | Main analysis: Baseline characteristics, T-learner, DUI, channel decomposition, behavioral and decomposition 2x2 |
| `sensitivity_analysis.R` | Total DUI re-estimated with four alternative base learners |
| `subgroup_analysis.R` | Primary care vs specialty care subgroup re-fit |

### Scope of the code released here

These scripts cover the **statistical analysis** only. The upstream
cohort-construction code, which queries the institutional electronic health
record and applies protected-health-information handling specific to our data
environment, is not included because it is inseparable from that environment and
cannot be executed outside it. The analysis scripts begin from a patient-level
analytic file described under **Expected input** below.

## Data availability

The individual-level electronic health record data underlying these analyses
were obtained from University of Pennsylvania and contain
protected health information. They cannot be publicly deposited. Requests for
access should be directed to yiwenlu@sas.upenn.edu and ychen123@pennmedicine.upenn.edu, subject to institutional approval and a
data use agreement.

## Expected input

Each script reads from a single directory set by `data_dir` at the top of the
file, and writes tables and figures to `results_dir`. Both are placeholders and
must be edited before running.

| Script | Reads | Writes |
|---|---|---|
| `analysis.R` | `telemedicine_cohort_enriched_<horizon>.csv` | `main_results_<horizon>.xlsx`, diagnostic PNGs, `telemedicine_analysis_final_<horizon>.csv` |
| `sensitivity_analysis.R` | `telemedicine_analysis_final_{3m,6m}.csv` | `sensitivity_models.xlsx` |
| `subgroup_analysis.R` | `telemedicine_analysis_final_{3m,6m}.csv`, provider-specialty crosswalk | `subgroup_results_{3m,6m}.xlsx` |

`analysis.R` must be run before the other two, because it produces
`telemedicine_analysis_final_<horizon>.csv`.

`subgroup_analysis.R` additionally requires a provider-specialty crosswalk,
supplied at the path set by `xwalk_path`. It is a CSV with two columns,
`specialty` and `care_category`, where `care_category` takes the values
`Primary Care` or `Specialty Care`. Index visits whose `PROV_SPECIALTY` does not
match a crosswalk row are reported in the `Unmatched_Specialty` sheet and
excluded from both subgroups.

Required columns in `telemedicine_cohort_enriched_<horizon>.csv`:

| Column | Description |
|---|---|
| `ENCOUNTER_TYPE` | Index visit modality label; telemedicine is defined by `tele_list` |
| `INDEX_CCSR` | CCSR category of the index visit primary diagnosis |
| `AGE_ENC`, `GENDER` | Age at index encounter, sex |
| `race_ethnicity` | Collapsed race and ethnicity category |
| `fin_class_category` | Financial class category |
| `CCI` | Charlson comorbidity index over the baseline lookback window |
| `n_inpatient_prior`, `n_ed_prior` | Baseline inpatient and emergency visit counts |
| `distance` | Distance between patient ZIP and place-of-service ZIP |
| `cohort_entry_month` | Index month, `YYYY-MM` |
| `n_followup_same_ccsr` | Same-CCSR follow-up visits within the horizon |
| `n_followup_tele`, `n_followup_f2f` | Same-CCSR follow-up visits by modality |
| `prior_same_ccsr_has_tele` | Baseline same-CCSR telemedicine use, for the behavioral 2x2 washout |
| `first_followup_modality` | Modality of the first same-CCSR follow-up visit |
| `PROV_SPECIALTY` | Index visit provider specialty, used by `subgroup_analysis.R` |

## Analysis summary

- **Cohort**: adult outpatient index visits, entry window 2024-01-01 to
  2024-11-30, restricted to the four high-virtualizability CCSR domains (END,
  MBD, MUS, SKN), complete cases on all covariates.
- **Diagnosis grouping**: ICD-10-CM diagnoses were mapped to CCSR categories
  during cohort construction using the AHRQ HCUP Clinical Classifications
  Software Refined (CCSR) for ICD-10-CM Diagnoses, reference file version
  2023-1, available from https://hcup-us.ahrq.gov. Each diagnosis was assigned
  to its first listed CCSR category. `CCSR_broad` is the three-letter body
  system prefix of `INDEX_CCSR`.
- **Horizons**: 3-month and 6-month follow-up on the identical cohort. Set
  `horizon` at the top of `analysis.R` and run once per horizon.
- **Estimator**: T-learner with XGBoost base learners (Poisson count objective),
  fit separately in the telemedicine and face-to-face arms on an 80/20 split
  (seed 42). DUI is the difference in predicted same-CCSR follow-up visits under
  the two counterfactual index modalities.
- **Channel decomposition**: four additional arm-specific models predict
  telemedicine-channel and face-to-face-channel follow-up separately, giving
  `eu_tele_channel` and `eu_f2f_channel`. An additivity gap against total DUI is
  reported.
- **Decomposition 2x2**: patients are assigned to one of four quadrants by
  thresholding both channel effects at delta. `analysis.R` uses delta = 0.25.
  `subgroup_analysis.R` uses delta = 0.25 at 3 months and a data-scaled delta at
  6 months, computed from the full-cohort ratio of mean total DUI across the two
  horizons.
- **Sensitivity**: total DUI re-estimated with XGBoost, Poisson GLM, random
  forest (ranger), and LightGBM base learners.
- **Subgroup**: the full T-learner is re-fit within primary care and specialty
  care separately, with a minimum of 200 patients per treatment arm.

## Requirements

R (>= 4.2) with:

```r
install.packages(c("data.table", "dplyr", "xgboost", "ggplot2",
                   "tableone", "openxlsx", "ranger", "lightgbm"))
```

`lightgbm` may require cmake on some systems.

## License

Released under the MIT License. See [LICENSE](LICENSE) for the full text.

The license covers the code in this repository only. It does not extend to the
underlying electronic health record data, which are not released, or to the AHRQ
HCUP CCSR reference file, which is distributed by AHRQ under its own terms of
use.
