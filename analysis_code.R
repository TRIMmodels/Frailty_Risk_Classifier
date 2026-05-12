# ============================================================
# TRIM frailty full analysis script
# Peng et al. Integrated Proteomic and Metabolomic Profiling
# for Frailty Identification and Prediction
#
# Purpose:
# This script provides a sequential workflow
# for data preparation, frailty index construction, omics
# preprocessing, feature selection, TRIM-FI/TRIM-FP modelling,
# external validation, model performance evaluation, and sensitivity
# analyses.
#
# Notes:
# - UK Biobank is used as the derivation cohort.
# - ESTHER is used as the external validation cohort.
# - Main TRIM models use 17 proteins + 9 metabolites selected
#   by bootstrap-enhanced LASSO for prevalent frailty in UKB.
# - Incident frailty is modelled using interval-censored
#   competing-risks regression with death as the competing event.
# ============================================================



# ============================================================
# 00. Project setup
# ============================================================

# ---- 00.1 Clear workspace ----
rm(list = ls())

# ---- 00.2 Set seed ----
set.seed(20260430)

# ---- 00.3 Required packages ----
required_packages <- c(
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "readr",
  "readxl",
  "writexl",
  "ggplot2",
  "glmnet",
  "missRanger",
  "pROC",
  "timeROC",
  "nricens",
  "rms",
  "rmda",
  "intccr",
  "survival",
  "broom",
  "broom.helpers",
  "scales"
)

new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

invisible(lapply(required_packages, library, character.only = TRUE))

# ---- 00.4 Project paths ----
# Modify root_dir to your local project folder.
root_dir <- getwd()

dir_data_raw      <- file.path(root_dir, "data_raw")
dir_data_processed <- file.path(root_dir, "data_processed")
dir_outputs       <- file.path(root_dir, "outputs")
dir_tables        <- file.path(dir_outputs, "tables")
dir_figures       <- file.path(dir_outputs, "figures")
dir_models        <- file.path(dir_outputs, "model_objects")
dir_logs          <- file.path(dir_outputs, "logs")

dirs_to_create <- c(
  dir_data_raw,
  dir_data_processed,
  dir_outputs,
  dir_tables,
  dir_figures,
  dir_models,
  dir_logs
)

invisible(lapply(dirs_to_create, dir.create, recursive = TRUE, showWarnings = FALSE))

# ---- 00.5 Global analysis settings ----
analysis_settings <- list(
  derivation_cohort = "UKB",
  validation_cohort = "ESTHER",

  # Missing data imputation, manuscript Methods:
  # random forest-based imputation with 100 trees and 5 iterations.
  imputation_num_trees = 100,
  imputation_max_iter = 5,
  imputation_pmm_k = 5,

  # Bootstrap-enhanced LASSO, manuscript Methods:
  # 1000 bootstrap samples, 10-fold CV, lambda.1se rule.
  lasso_bootstrap_B = 1000,
  lasso_nfolds = 10,
  lasso_lambda_rule = "lambda.1se",
  lasso_selection_threshold = 0.95,

  # Multiple testing correction.
  p_adjust_method = "BH",

  # Risk categories used for categorical NRI:
  # <5%, 5–10%, >10%.
  nri_cutoffs = c(0.05, 0.10)
)

# ---- 00.6 Helper: display names for output ----
make_display_name <- function(x) {
  dplyr::recode(
    x,
    "LA_pct" = "LA-pct",
    "ApoB_by_ApoA1" = "ApoB-by-ApoA1",
    "XL_HDL_FC" = "XL-HDL-FC",
    "S_HDL_PL" = "S-HDL-PL",
    .default = x
  )
}


# ============================================================
# 01. Biomarker panels and covariates
# ============================================================

# ---- 01.1 Main analysis biomarker panel ----
# Main TRIM models use 17 proteins + 9 metabolites selected from
# prevalent frailty LASSO in UKB.

main_analysis <- list(
  proteins = c(
    "TNFRSF11B", "DNER", "FGF23", "FGF21", "HGF", "VEGFA",
    "CXCL8", "CXCL10", "CCL11", "CCL20", "CDCP1", "CD5",
    "TGFA", "TGFB1", "TNFRSF9", "IL18", "IL6"
  ),

  metabolites = c(
    "His", "Glucose", "Lactate", "Albumin", "GlycA",
    "Gln", "LA_pct", "ApoB_by_ApoA1", "XL_HDL_FC"
  )
)

main_analysis$omics <- c(
  main_analysis$proteins,
  main_analysis$metabolites
)

# ---- 01.2 Sensitivity analysis biomarker panel ----
# Incident-frailty-derived panel from Supplementary Table 7.
# Note:
# Supplementary Table 7 text says "6 metabolites", but the footnote
# and Supplementary Table 8 list 5 metabolites.

sensitivity_analysis <- list(
  proteins = c(
    "CCL3", "FGF23", "HGF", "IL12B", "IL6",
    "LIFR", "TGFA", "TNFRSF9", "VEGFA"
  ),

  metabolites = c(
    "Albumin", "Glucose", "GlycA", "Lactate", "S_HDL_PL"
  )
)

sensitivity_analysis$omics <- c(
  sensitivity_analysis$proteins,
  sensitivity_analysis$metabolites
)

# ---- 01.3 Conventional risk factors ----
# Reference model covariates:
# age, sex, education, smoking status, physical activity,
# alcohol consumption, BMI, dyslipidemia, hypertension, diabetes,
# cardiovascular disease, and history of cancer.

continuous_covariates <- c("age")

categorical_covariates <- c(
  "sex",
  "education",
  "smoking",
  "physical_activity",
  "alcohol",
  "BMI_cat",
  "dyslipidemia",
  "hypertension",
  "diabetes",
  "CVD",
  "cancer"
)

reference_covariates <- c(
  continuous_covariates,
  categorical_covariates
)

# ---- 01.4 Expected category coding ----
# These levels are used to ensure consistent reference categories
# across UKB and ESTHER.

category_levels <- list(
  sex = c("Women", "Men"),
  education = c("Medium or low", "High"),
  smoking = c("Never smoker", "Former smoker", "Current smoker"),
  physical_activity = c("Low", "Moderate", "High"),
  alcohol = c("Abstainer", "Low", "Medium or high"),
  BMI_cat = c("<21.5", "21.5-<25", "25-<30", "30-<35", ">=35"),
  dyslipidemia = c("No", "Yes"),
  hypertension = c("No", "Yes"),
  diabetes = c("No", "Yes"),
  CVD = c("No", "Yes"),
  cancer = c("No", "Yes")
)

# ---- 01.5 Biomarker display names ----
# Data columns use syntactic names.

biomarker_display <- c(
  LA_pct = "LA-pct",
  ApoB_by_ApoA1 = "ApoB-by-ApoA1",
  XL_HDL_FC = "XL-HDL-FC",
  S_HDL_PL = "S-HDL-PL"
)

# ---- 01.6 Compatibility objects ----
# These objects keep the script compatible with earlier modular scripts.

proteins_main <- main_analysis$proteins
metabolites_main <- main_analysis$metabolites
omics_main <- main_analysis$omics

proteins_incident_sensitivity <- sensitivity_analysis$proteins
metabolites_incident_sensitivity <- sensitivity_analysis$metabolites
omics_incident_sensitivity <- sensitivity_analysis$omics


# ============================================================
# 02. Data loading
# ============================================================

# ---- 02.1 Expected input datasets ----
# This script assumes that harmonised participant-level UKB and ESTHER
# datasets have already been prepared.
#
# Required baseline datasets:
#   1. UKB baseline dataset
#   2. ESTHER baseline dataset
#
# Required follow-up datasets:
#   3. UKB follow-up dataset
#   4. ESTHER follow-up dataset
#
# Each baseline dataset should contain:
#   - id
#   - conventional risk factors
#   - 31 harmonised frailty index items: FI_item_01 ... FI_item_31
#   - proteomic biomarkers
#   - metabolomic biomarkers
#
# Each follow-up dataset should contain:
#   - id
#   - follow-up frailty index or follow-up FI items
#   - follow-up visit time
#   - death indicator, if applicable

# ---- 02.2 File paths ----
# Replace filenames with your actual harmonised data files.

file_ukb_baseline <- file.path(dir_data_raw, "ukb_baseline_harmonised.csv")
file_ukb_followup <- file.path(dir_data_raw, "ukb_followup_harmonised.csv")

file_esther_baseline <- file.path(dir_data_raw, "esther_baseline_harmonised.csv")
file_esther_followup <- file.path(dir_data_raw, "esther_followup_harmonised.csv")

# ---- 02.3 Read data ----

ukb_bl <- readr::read_csv(file_ukb_baseline, show_col_types = FALSE)
ukb_fu <- readr::read_csv(file_ukb_followup, show_col_types = FALSE)

esther_bl <- readr::read_csv(file_esther_baseline, show_col_types = FALSE)
esther_fu <- readr::read_csv(file_esther_followup, show_col_types = FALSE)

# ---- 02.4 Required baseline variables ----

frailty_deficit_vars <- paste0("FI_item_", sprintf("%02d", 1:31))

required_baseline_vars <- c(
  "id",
  reference_covariates,
  frailty_deficit_vars,
  omics_main
)

check_required_vars <- function(data, vars, data_name = "dataset") {
  missing_vars <- setdiff(vars, names(data))

  if (length(missing_vars) > 0) {
    stop(
      paste0(
        "The following required variables are missing in ",
        data_name,
        ": ",
        paste(missing_vars, collapse = ", ")
      )
    )
  }

  invisible(TRUE)
}

check_required_vars(ukb_bl, required_baseline_vars, "UKB baseline")
check_required_vars(esther_bl, required_baseline_vars, "ESTHER baseline")

# ---- 02.5 Harmonise categorical variables ----
# Reference categories are fixed according to the manuscript models.

apply_category_levels <- function(data, category_levels) {
  for (v in names(category_levels)) {
    if (v %in% names(data)) {
      data[[v]] <- factor(data[[v]], levels = category_levels[[v]])
    }
  }

  data
}

ukb_bl <- apply_category_levels(ukb_bl, category_levels)
esther_bl <- apply_category_levels(esther_bl, category_levels)

# ---- 02.6 Add cohort labels ----

ukb_bl <- ukb_bl %>%
  mutate(cohort = "UKB")

esther_bl <- esther_bl %>%
  mutate(cohort = "ESTHER")


# ============================================================
# 03. Frailty index construction
# ============================================================

# ---- 03.1 Frailty index definition ----
# Frailty was assessed using the Rockwood deficit accumulation approach.
#
# The manuscript used the same 31 health deficits in UKB and ESTHER:
#   1-7   disease history
#   8-11  major disease events
#   12-21 medications
#   22-27 functional limitations / ADL difficulties
#   28    self-rated general health
#   29-31 lifestyle-related factors
#
# Each deficit is scored on a 0-to-1 scale, where 0 indicates absence
# of the deficit and 1 indicates full presence; intermediate values are
# used for ordinal or graded deficits where applicable.
# FI is calculated as the sum of the 31 deficit scores divided by 31.

# ---- 03.2 Frailty thresholds ----
# Cohort-specific validated FI thresholds used in the manuscript:
#   UKB:    frailty defined as FI >= 0.51
#   ESTHER: frailty defined as FI > 0.35

frailty_thresholds <- list(
  UKB = list(value = 0.51, operator = ">="),
  ESTHER = list(value = 0.35, operator = ">")
)

classify_frailty <- function(FI, cohort) {
  if (cohort == "UKB") {
    return(as.integer(FI >= frailty_thresholds$UKB$value))
  }

  if (cohort == "ESTHER") {
    return(as.integer(FI > frailty_thresholds$ESTHER$value))
  }

  stop("Unknown cohort: ", cohort)
}

# ---- 03.3 Check frailty item values ----

check_frailty_item_values <- function(data, data_name = "dataset") {
  item_range <- data %>%
    summarise(
      across(
        all_of(frailty_deficit_vars),
        list(
          min = ~ min(.x, na.rm = TRUE),
          max = ~ max(.x, na.rm = TRUE)
        )
      )
    )

  invalid_items <- names(data)[
    names(data) %in% frailty_deficit_vars &
      purrr::map_lgl(
        data[names(data) %in% frailty_deficit_vars],
        ~ any(!is.na(.x) & (.x < 0 | .x > 1))
      )
  ]

  if (length(invalid_items) > 0) {
    stop(
      paste0(
        "Frailty item values outside [0, 1] detected in ",
        data_name,
        ": ",
        paste(invalid_items, collapse = ", ")
      )
    )
  }

  invisible(item_range)
}

check_frailty_item_values(ukb_bl, "UKB baseline")
check_frailty_item_values(esther_bl, "ESTHER baseline")

# ---- 03.4 Construct frailty index ----
# Conservative public implementation:
# FI is calculated as FI_sum / 31 when all 31 FI items are available.
#
# To avoid treating missing deficit values as zero, FI is set to missing
# unless all 31 FI items are available.
#
# If the manuscript pipeline imputed or otherwise handled missing FI items
# before FI construction, apply that preprocessing before this step.

derive_frailty_index_31items <- function(data) {
  data %>%
    mutate(
      across(
        all_of(frailty_deficit_vars),
        ~ as.numeric(.x)
      )
    ) %>%
    mutate(
      FI_n_available = rowSums(!is.na(across(all_of(frailty_deficit_vars)))),
      FI_sum = rowSums(across(all_of(frailty_deficit_vars)), na.rm = TRUE),
      FI = ifelse(FI_n_available == 31, FI_sum / 31, NA_real_)
    )
}

ukb_bl <- derive_frailty_index_31items(ukb_bl)
esther_bl <- derive_frailty_index_31items(esther_bl)

# ---- 03.5 Define prevalent frailty ----

ukb_bl <- ukb_bl %>%
  mutate(
    frailty_prevalent = classify_frailty(FI, "UKB")
  )

esther_bl <- esther_bl %>%
  mutate(
    frailty_prevalent = classify_frailty(FI, "ESTHER")
  )

# ---- 03.6 Prevalent frailty check ----
# Expected manuscript numbers after applying all inclusion/exclusion criteria:
#   UKB:    n = 19,151; frail = 1,863; 9.7%
#   ESTHER: n = 5,031;  frail = 494;   9.8%

frailty_check <- bind_rows(
  ukb_bl %>%
    summarise(
      cohort = "UKB",
      n = n(),
      frail_n = sum(frailty_prevalent == 1, na.rm = TRUE),
      frail_percent = 100 * mean(frailty_prevalent == 1, na.rm = TRUE),
      FI_mean = mean(FI, na.rm = TRUE),
      FI_sd = sd(FI, na.rm = TRUE),
      FI_missing_n = sum(is.na(FI))
    ),

  esther_bl %>%
    summarise(
      cohort = "ESTHER",
      n = n(),
      frail_n = sum(frailty_prevalent == 1, na.rm = TRUE),
      frail_percent = 100 * mean(frailty_prevalent == 1, na.rm = TRUE),
      FI_mean = mean(FI, na.rm = TRUE),
      FI_sd = sd(FI, na.rm = TRUE),
      FI_missing_n = sum(is.na(FI))
    )
)

print(frailty_check)

readr::write_csv(
  frailty_check,
  file.path(dir_logs, "frailty_prevalence_check.csv")
)

# ---- 03.7 Optional recoding template for FI items ----
# Use this section only if your datasets do not already contain
# FI_item_01 ... FI_item_31.
#
# Important:
# The exact raw variable names differ between UKB and ESTHER.
# Therefore, raw variables should first be harmonised into the
# 31 FI items below.

# Example template:
#
# data <- data %>%
#   mutate(
#     FI_item_01 = coronary_artery_disease,
#     FI_item_02 = heart_failure,
#     FI_item_03 = diabetes_history,
#     FI_item_04 = cancer_history,
#     FI_item_05 = glaucoma,
#     FI_item_06 = cataract,
#     FI_item_07 = parkinsons_disease,
#     FI_item_08 = myocardial_infarction,
#     FI_item_09 = stroke,
#     FI_item_10 = joint_replacement,
#     FI_item_11 = hip_fracture,
#     FI_item_12 = antihypertensive_drugs,
#     FI_item_13 = lipid_lowering_drugs,
#     FI_item_14 = vasodilators,
#     FI_item_15 = cardiac_glycosides,
#     FI_item_16 = prescribed_aspirin,
#     FI_item_17 = anti_osteoporotic_drugs,
#     FI_item_18 = anxiolytics,
#     FI_item_19 = sedatives,
#     FI_item_20 = anti_dementia_drugs,
#     FI_item_21 = prostate_hyperplasia_or_incontinence_drugs,
#     FI_item_22 = difficulty_moderate_activities,
#     FI_item_23 = difficulty_climbing_stairs,
#     FI_item_24 = pain_related_limitation,
#     FI_item_25 = physical_health_limits_activities,
#     FI_item_26 = limited_social_contacts,
#     FI_item_27 = emotional_problem_accomplished_less,
#     FI_item_28 = case_when(
#       self_rated_health == "Excellent" ~ 0,
#       self_rated_health == "Good" ~ 0.333,
#       self_rated_health == "Fair" ~ 0.667,
#       self_rated_health == "Poor" ~ 1,
#       TRUE ~ NA_real_
#     ),
#     FI_item_29 = ifelse(BMI < 20, 1, 0),
#     FI_item_30 = case_when(
#       BMI >= 35 ~ 1,
#       BMI >= 30 & BMI < 35 ~ 0.5,
#       BMI < 30 ~ 0,
#       TRUE ~ NA_real_
#     ),
#     FI_item_31 = ifelse(vigorous_physical_activity_hours_per_week == 0, 1, 0)
#   )


# ============================================================
# 04. Omics preprocessing
# ============================================================

# ---- 04.1 Omics preprocessing rules ----
#
# Proteomics:
#   - Protein abundance values are NPX values on log2 scale.
#   - Proteins are standardised to z-scores within each cohort.
#
# Metabolomics:
#   - Metabolite concentrations are log2-transformed.
#   - Then standardised to z-scores within each cohort.
#
# Important:
#   - Standardisation is performed separately in UKB and ESTHER.
#   - This avoids using information from the validation cohort to
#     scale the derivation cohort, and matches the manuscript approach.
#   - Main analyses use 67 common inflammation-related proteins and
#     249 metabolites after QC, but the TRIM models use the selected
#     17 proteins + 9 metabolites defined in Section 01.

# ---- 04.2 Define variables to preprocess ----

proteomics_vars_to_scale <- unique(c(
  main_analysis$proteins,
  sensitivity_analysis$proteins
))

metabolomics_vars_to_log_scale <- unique(c(
  main_analysis$metabolites,
  sensitivity_analysis$metabolites
))

omics_vars_to_check <- unique(c(
  proteomics_vars_to_scale,
  metabolomics_vars_to_log_scale
))

# ---- 04.3 Check availability of omics variables ----

check_omics_vars <- function(data, vars, data_name = "dataset") {
  missing_vars <- setdiff(vars, names(data))

  if (length(missing_vars) > 0) {
    warning(
      paste0(
        "The following omics variables are missing in ",
        data_name,
        ": ",
        paste(missing_vars, collapse = ", ")
      )
    )
  }

  invisible(intersect(vars, names(data)))
}

ukb_omics_available <- check_omics_vars(
  ukb_bl,
  omics_vars_to_check,
  "UKB baseline"
)

esther_omics_available <- check_omics_vars(
  esther_bl,
  omics_vars_to_check,
  "ESTHER baseline"
)

# ---- 04.4 Helper: safe log2 transformation ----
# Metabolite values should be positive before log2 transformation.
# If zero or negative values exist, stop and check the raw data.
# The manuscript used log2-transformed metabolomic biomarker concentrations.

safe_log2 <- function(x, var_name = "variable") {
  x_num <- suppressWarnings(as.numeric(x))

  if (any(!is.na(x_num) & x_num <= 0)) {
    stop(
      paste0(
        "Non-positive values detected in ",
        var_name,
        ". Log2 transformation requires values > 0."
      )
    )
  }

  log2(x_num)
}

# ---- 04.5 Helper: z-score standardisation ----

zscore <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  as.numeric(scale(x_num))
}

# ---- 04.6 Preprocess omics within one cohort ----

preprocess_omics_one_cohort <- function(
    data,
    proteins,
    metabolites,
    cohort_name = "cohort"
) {
  proteins_present <- intersect(proteins, names(data))
  metabolites_present <- intersect(metabolites, names(data))

  data_out <- data

  # Proteomics:
  # NPX values are already on log2 scale, so only z-score standardisation.
  if (length(proteins_present) > 0) {
    data_out <- data_out %>%
      mutate(
        across(
          all_of(proteins_present),
          ~ zscore(.x)
        )
      )
  }

  # Metabolomics:
  # log2 transformation followed by z-score standardisation.
  if (length(metabolites_present) > 0) {
    for (v in metabolites_present) {
      data_out[[v]] <- safe_log2(data_out[[v]], var_name = paste0(cohort_name, "::", v))
      data_out[[v]] <- zscore(data_out[[v]])
    }
  }

  data_out
}

# ---- 04.7 Apply omics preprocessing separately in UKB and ESTHER ----

ukb_bl <- preprocess_omics_one_cohort(
  data = ukb_bl,
  proteins = proteomics_vars_to_scale,
  metabolites = metabolomics_vars_to_log_scale,
  cohort_name = "UKB baseline"
)

esther_bl <- preprocess_omics_one_cohort(
  data = esther_bl,
  proteins = proteomics_vars_to_scale,
  metabolites = metabolomics_vars_to_log_scale,
  cohort_name = "ESTHER baseline"
)

# ---- 04.8 Save z-score summaries for checking ----

summarise_omics_zscores <- function(data, vars, cohort_name) {
  vars_present <- intersect(vars, names(data))

  data %>%
    summarise(
      across(
        all_of(vars_present),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          missing = ~ sum(is.na(.x))
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    mutate(cohort = cohort_name, .before = 1)
}

omics_zscore_check <- bind_rows(
  summarise_omics_zscores(ukb_bl, omics_vars_to_check, "UKB"),
  summarise_omics_zscores(esther_bl, omics_vars_to_check, "ESTHER")
)

readr::write_csv(
  omics_zscore_check,
  file.path(dir_logs, "omics_zscore_check.csv")
)


# ============================================================
# 05. Missing data imputation
# ============================================================

# ---- 05.1 Imputation rules ----
#
# Missing values were imputed using a random forest-based imputation method,
# performed separately within each cohort.
#
# Parameters:
#   - num.trees = 100
#   - maxiter   = 5
#
# R package:
#   - missRanger
#
# Important:
#   - UKB and ESTHER are imputed separately.
#   - Baseline and follow-up analysis datasets should also be handled
#     separately if their available variables differ.
#   - The outcome variables should not be imputed.
#   - Participant ID and cohort labels should not be imputed.

# ---- 05.2 Variables used for baseline imputation ----
# For baseline analyses, impute covariates and selected omics biomarkers.
# The prevalent frailty outcome is excluded from the imputation target.

baseline_imputation_vars <- unique(c(
  reference_covariates,
  omics_main,
  "frailty_prevalent"
))

baseline_imputation_vars_ukb <- intersect(
  baseline_imputation_vars,
  names(ukb_bl)
)

baseline_imputation_vars_esther <- intersect(
  baseline_imputation_vars,
  names(esther_bl)
)

# ---- 05.3 Helper: prepare imputation dataset ----

prepare_imputation_data <- function(
    data,
    vars,
    outcome_vars = c("frailty_prevalent"),
    id_vars = c("id", "cohort")
) {
  vars_present <- intersect(vars, names(data))

  impute_data <- data %>%
    select(all_of(c(intersect(id_vars, names(data)), vars_present)))

  # Convert character variables to factors for missRanger.
  impute_data <- impute_data %>%
    mutate(
      across(
        where(is.character),
        as.factor
      )
    )

  # Ensure outcome variables are factors or numeric but retained only as predictors.
  # They will be restored after imputation to avoid changing the observed outcome.
  impute_data
}

# ---- 05.4 Helper: impute one cohort ----

impute_one_cohort <- function(
    data,
    imputation_vars,
    outcome_vars = c("frailty_prevalent"),
    id_vars = c("id", "cohort"),
    num_trees = analysis_settings$imputation_num_trees,
    max_iter = analysis_settings$imputation_max_iter,
    pmm_k = analysis_settings$imputation_pmm_k,
    seed = 20260430
) {
  set.seed(seed)

  impute_data <- prepare_imputation_data(
    data = data,
    vars = imputation_vars,
    outcome_vars = outcome_vars,
    id_vars = id_vars
  )

  # Store non-imputed identifiers and outcomes.
  id_data <- impute_data %>%
    select(any_of(id_vars))

  outcome_data <- impute_data %>%
    select(any_of(outcome_vars))

  # Variables eligible for imputation:
  # exclude ID variables and outcomes.
  impute_targets <- setdiff(
    names(impute_data),
    c(id_vars, outcome_vars)
  )

  dat_for_imputation <- impute_data %>%
    select(all_of(impute_targets))

  # missRanger imputes all missing values in dat_for_imputation.
  #
  # Predictive mean matching is used with pmm.k = 5 to preserve
  # the empirical distribution of imputed values. If the original
  # analysis pipeline used a different pmm.k setting, modify pmm_k
  # accordingly.
  dat_imputed <- missRanger::missRanger(
    data = dat_for_imputation,
    num.trees = num_trees,
    maxiter = max_iter,
    pmm.k = pmm_k,
    seed = seed,
    verbose = 1
  )

  out <- bind_cols(
    id_data,
    dat_imputed,
    outcome_data
  )

  # Restore variables not included in imputation.
  untouched_vars <- setdiff(names(data), names(out))

  out <- data %>%
    select(any_of(untouched_vars)) %>%
    bind_cols(out)

  # Re-apply categorical levels after imputation.
  out <- apply_category_levels(out, category_levels)

  out
}

# ---- 05.5 Baseline imputation in UKB and ESTHER ----

ukb_bl_imp <- impute_one_cohort(
  data = ukb_bl,
  imputation_vars = baseline_imputation_vars_ukb,
  outcome_vars = c("frailty_prevalent"),
  id_vars = c("id", "cohort"),
  seed = 20260430
)

esther_bl_imp <- impute_one_cohort(
  data = esther_bl,
  imputation_vars = baseline_imputation_vars_esther,
  outcome_vars = c("frailty_prevalent"),
  id_vars = c("id", "cohort"),
  seed = 20260430
)

# ---- 05.6 Check missingness before and after imputation ----

missingness_summary <- function(data_before, data_after, vars, cohort_name) {
  vars_present <- intersect(vars, names(data_before))
  vars_present <- intersect(vars_present, names(data_after))

  if (length(vars_present) == 0) {
    return(
      tibble(
        cohort = character(),
        variable = character(),
        missing_before = integer(),
        missing_after = integer(),
        missing_percent_before = numeric(),
        missing_percent_after = numeric()
      )
    )
  }

  tibble(
    cohort = cohort_name,
    variable = vars_present,
    missing_before = purrr::map_int(
      vars_present,
      ~ sum(is.na(data_before[[.x]]))
    ),
    missing_after = purrr::map_int(
      vars_present,
      ~ sum(is.na(data_after[[.x]]))
    ),
    missing_percent_before = purrr::map_dbl(
      vars_present,
      ~ 100 * mean(is.na(data_before[[.x]]))
    ),
    missing_percent_after = purrr::map_dbl(
      vars_present,
      ~ 100 * mean(is.na(data_after[[.x]]))
    )
  )
}

baseline_missingness_check <- bind_rows(
  missingness_summary(
    data_before = ukb_bl,
    data_after = ukb_bl_imp,
    vars = baseline_imputation_vars_ukb,
    cohort_name = "UKB"
  ),
  missingness_summary(
    data_before = esther_bl,
    data_after = esther_bl_imp,
    vars = baseline_imputation_vars_esther,
    cohort_name = "ESTHER"
  )
)

readr::write_csv(
  baseline_missingness_check,
  file.path(dir_logs, "baseline_missingness_imputation_check.csv")
)

# ---- 05.7 Save processed baseline datasets ----

readr::write_csv(
  ukb_bl_imp,
  file.path(dir_data_processed, "ukb_baseline_imputed_z.csv")
)

readr::write_csv(
  esther_bl_imp,
  file.path(dir_data_processed, "esther_baseline_imputed_z.csv")
)

# ---- 05.8 Notes for follow-up imputation ----
# Follow-up analysis datasets will be constructed in Section 08.
# If follow-up-specific variables have missing values, the same imputation
# principle should be applied:
#
#   - UKB and ESTHER separately
#   - num.trees = 100
#   - maxiter = 5
#   - do not impute incident frailty outcome, death indicator, or time intervals
#
# This avoids leakage from outcome/event information into predictors.



# ============================================================
# 06. Main analysis: prevalent frailty LASSO
# ============================================================

# ---- 06.1 Analysis objective ----
# Identify a parsimonious set of proteomic and metabolomic biomarkers
# associated with prevalent frailty at baseline in UKB.
#
#  LASSO regression:
#   - Outcome: prevalent frailty at baseline
#   - Cohort: UKB derivation cohort
#   - Biomarkers: proteomic and metabolomic biomarkers jointly entered
#   - Adjustment: reference covariates always included
#   - Bootstrap samples: 1000
#   - Cross-validation: 10-fold
#   - Lambda rule: lambda.1se
#   - Selection threshold: selected in >=95% bootstrap samples
#
# Important:
# In the public/demo script, omics_main already contains the final
# LASSO-selected 17 proteins + 9 metabolites.
# If full candidate omics data are available, replace lasso_candidate_omics
# with the full QC-passed candidate set:
#   - 67 inflammation-related proteins
#   - 249 metabolites

# ---- 06.2 Candidate biomarkers for LASSO ----
# For exact reproduction from raw candidate features, use all QC-passed
# candidate biomarkers. Here, omics_main is used as the final panel
# selected panel for demonstration and downstream consistency.

# Public/example implementation:
# The LASSO-selected 26 biomarkers are used here for downstream
# model fitting and demonstration.
#
# To reproduce the original discovery LASSO exactly, replace
# lasso_candidate_omics with all QC-passed candidate biomarkers:
#   - 67 inflammation-related proteins
#   - 249 metabolites

lasso_candidate_omics <- omics_main

# ---- 06.3 Prepare model matrix for LASSO ----
# Reference covariates are forced into the model by setting penalty.factor = 0.
# Omics biomarkers are penalised by setting penalty.factor = 1.

make_lasso_design <- function(
    data,
    outcome,
    forced_covariates,
    penalised_biomarkers
) {
  vars_needed <- unique(c(outcome, forced_covariates, penalised_biomarkers))

  dat <- data %>%
    select(all_of(vars_needed)) %>%
    filter(!is.na(.data[[outcome]]))

  # Ensure categorical variables are factors with manuscript reference levels.
  dat <- apply_category_levels(dat, category_levels)

  formula_x <- as.formula(
    paste(
      "~",
      paste(c(forced_covariates, penalised_biomarkers), collapse = " + ")
    )
  )

  x <- model.matrix(formula_x, data = dat)[, -1, drop = FALSE]
  y <- as.integer(dat[[outcome]])

  if (!all(y %in% c(0, 1))) {
    stop("LASSO outcome must be coded as 0/1.")
  }

  # Penalise only biomarker columns.
  # Because categorical covariates expand into dummy variables, forced
  # covariate columns are identified by their variable-name prefixes.
  forced_pattern <- paste0("^(", paste(forced_covariates, collapse = "|"), ")")
  biomarker_pattern <- paste0("^(", paste(penalised_biomarkers, collapse = "|"), ")$")

  penalty_factor <- rep(1, ncol(x))
  penalty_factor[grepl(forced_pattern, colnames(x))] <- 0
  penalty_factor[grepl(biomarker_pattern, colnames(x))] <- 1

  list(
    data = dat,
    x = x,
    y = y,
    penalty_factor = penalty_factor,
    colnames = colnames(x)
  )
}

lasso_design_ukb <- make_lasso_design(
  data = ukb_bl_imp,
  outcome = "frailty_prevalent",
  forced_covariates = reference_covariates,
  penalised_biomarkers = lasso_candidate_omics
)

# ---- 06.4 Bootstrap-enhanced LASSO function ----

run_bootstrap_lasso <- function(
    x,
    y,
    penalty_factor,
    penalised_biomarkers,
    B = analysis_settings$lasso_bootstrap_B,
    nfolds = analysis_settings$lasso_nfolds,
    lambda_rule = analysis_settings$lasso_lambda_rule,
    selection_threshold = analysis_settings$lasso_selection_threshold,
    seed = 20260430
) {
  set.seed(seed)

  n <- nrow(x)

  selection_matrix <- matrix(
    0,
    nrow = B,
    ncol = length(penalised_biomarkers)
  )

  colnames(selection_matrix) <- penalised_biomarkers

  lambda_values <- numeric(B)

  for (b in seq_len(B)) {
    if (b %% 50 == 0) {
      message("Bootstrap LASSO iteration: ", b, " / ", B)
    }

    idx <- sample(seq_len(n), size = n, replace = TRUE)

    x_b <- x[idx, , drop = FALSE]
    y_b <- y[idx]

    cvfit <- glmnet::cv.glmnet(
      x = x_b,
      y = y_b,
      family = "binomial",
      alpha = 1,
      nfolds = nfolds,
      penalty.factor = penalty_factor,
      standardize = FALSE,
      type.measure = "deviance"
    )

    lambda_use <- if (lambda_rule == "lambda.1se") {
      cvfit$lambda.1se
    } else if (lambda_rule == "lambda.min") {
      cvfit$lambda.min
    } else {
      stop("Unknown lambda_rule: ", lambda_rule)
    }

    lambda_values[b] <- lambda_use

    coef_b <- as.matrix(coef(cvfit, s = lambda_use))
    selected_terms <- rownames(coef_b)[as.numeric(coef_b[, 1]) != 0]
    selected_terms <- setdiff(selected_terms, "(Intercept)")

    # Count biomarker selection only.
    for (v in penalised_biomarkers) {
      selection_matrix[b, v] <- as.integer(v %in% selected_terms)
    }
  }

  selection_frequency <- colMeans(selection_matrix)

  selection_summary <- tibble(
    biomarker = names(selection_frequency),
    selection_frequency = as.numeric(selection_frequency),
    selected = selection_frequency >= selection_threshold
  ) %>%
    arrange(desc(selection_frequency), biomarker)

  selected_biomarkers <- selection_summary %>%
    filter(selected) %>%
    pull(biomarker)

  list(
    selection_matrix = selection_matrix,
    selection_summary = selection_summary,
    selected_biomarkers = selected_biomarkers,
    lambda_values = lambda_values
  )
}

# ---- 06.5 Run or load bootstrap LASSO ----
# Running B = 1000 can take time. If a saved result exists, load it.

file_lasso_prevalent <- file.path(
  dir_models,
  "ukb_prevalent_bootstrap_lasso.rds"
)

if (file.exists(file_lasso_prevalent)) {
  lasso_prevalent_ukb <- readRDS(file_lasso_prevalent)
} else {
  lasso_prevalent_ukb <- run_bootstrap_lasso(
    x = lasso_design_ukb$x,
    y = lasso_design_ukb$y,
    penalty_factor = lasso_design_ukb$penalty_factor,
    penalised_biomarkers = lasso_candidate_omics,
    B = analysis_settings$lasso_bootstrap_B,
    nfolds = analysis_settings$lasso_nfolds,
    lambda_rule = analysis_settings$lasso_lambda_rule,
    selection_threshold = analysis_settings$lasso_selection_threshold,
    seed = 20260430
  )

  saveRDS(lasso_prevalent_ukb, file_lasso_prevalent)
}

# ---- 06.6 Save LASSO selection summary ----

readr::write_csv(
  lasso_prevalent_ukb$selection_summary,
  file.path(dir_tables, "ukb_prevalent_lasso_selection_summary.csv")
)

# ---- 06.7 LASSO-selected biomarkers panel----
# The final main TRIM panel is fixed to the LASSO-selected
# 17 proteins + 9 metabolites from Section 01.

selected_omics_main <- omics_main

# Optional consistency check:
# If full candidate LASSO was run, selected_biomarkers should match
# selected_omics_main after harmonising names.

lasso_selected_from_run <- lasso_prevalent_ukb$selected_biomarkers

if (length(lasso_selected_from_run) > 0) {
  missing_from_lasso <- setdiff(selected_omics_main, lasso_selected_from_run)

  if (length(missing_from_lasso) > 0) {
    warning(
      paste0(
        "Some biomarkers were not selected in the current ",
        "LASSO run. This is expected if lasso_candidate_omics was set to ",
        "the final panel for demonstration or if data differ from the ",
        "manuscript analysis: ",
        paste(missing_from_lasso, collapse = ", ")
      )
    )
  }
}


# ============================================================
# 07. TRIM-FI: baseline frailty identification
# ============================================================

# ---- 07.1 Analysis objective ----
# Fit logistic regression models for prevalent frailty at baseline.
#
# Models:
#   1. Reference model:
#      conventional risk factors only
#
#   2. Reference + Proteomics:
#      conventional risk factors + 17 selected proteins
#
#   3. Reference + Metabolomics:
#      conventional risk factors + 9 selected metabolites
#
#   4. TRIM-FI:
#      conventional risk factors + 17 proteins + 9 metabolites
#
# Cohort:
#   - UKB is the derivation cohort.
#   - ESTHER external validation is performed in Section 09 by applying
#     UKB-derived coefficients directly, without refitting the main model.

# ---- 07.2 Define TRIM-FI model variable sets ----

trim_fi_variable_sets <- list(
  Reference = reference_covariates,
  Reference_plus_Proteomics = c(reference_covariates, main_analysis$proteins),
  Reference_plus_Metabolomics = c(reference_covariates, main_analysis$metabolites),
  TRIM_FI = c(reference_covariates, main_analysis$omics)
)

# ---- 07.3 Helper: build logistic formula ----

make_model_formula <- function(outcome, predictors) {
  as.formula(
    paste(
      outcome,
      "~",
      paste(predictors, collapse = " + ")
    )
  )
}

# ---- 07.4 Helper: fit logistic regression ----

fit_logistic_model <- function(data, outcome, predictors) {
  vars_needed <- unique(c(outcome, predictors))

  dat <- data %>%
    select(all_of(vars_needed)) %>%
    filter(!is.na(.data[[outcome]]))

  dat <- apply_category_levels(dat, category_levels)

  formula <- make_model_formula(outcome, predictors)

  glm(
    formula = formula,
    data = dat,
    family = binomial(link = "logit")
  )
}

# ---- 07.5 Fit TRIM-FI models in UKB ----

trim_fi_models_ukb <- list(
  Reference = fit_logistic_model(
    data = ukb_bl_imp,
    outcome = "frailty_prevalent",
    predictors = trim_fi_variable_sets$Reference
  ),

  Reference_plus_Proteomics = fit_logistic_model(
    data = ukb_bl_imp,
    outcome = "frailty_prevalent",
    predictors = trim_fi_variable_sets$Reference_plus_Proteomics
  ),

  Reference_plus_Metabolomics = fit_logistic_model(
    data = ukb_bl_imp,
    outcome = "frailty_prevalent",
    predictors = trim_fi_variable_sets$Reference_plus_Metabolomics
  ),

  TRIM_FI = fit_logistic_model(
    data = ukb_bl_imp,
    outcome = "frailty_prevalent",
    predictors = trim_fi_variable_sets$TRIM_FI
  )
)

saveRDS(
  trim_fi_models_ukb,
  file.path(dir_models, "trim_fi_models_ukb.rds")
)

# ---- 07.6 Extract TRIM-FI coefficients ----

extract_glm_coefficients <- function(model, model_name) {
  broom::tidy(model) %>%
    mutate(
      model = model_name,
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std.error),
      CI_high = exp(estimate + 1.96 * std.error),
      p_value = p.value
    ) %>%
    select(
      model,
      term,
      estimate,
      std.error,
      OR,
      CI_low,
      CI_high,
      p_value
    )
}

trim_fi_coefficients_ukb <- purrr::imap_dfr(
  trim_fi_models_ukb,
  extract_glm_coefficients
)

readr::write_csv(
  trim_fi_coefficients_ukb,
  file.path(dir_tables, "trim_fi_coefficients_ukb.csv")
)

# ---- 07.7 Single-biomarker prevalent frailty associations ----
# Each selected biomarker is modelled separately, adjusted for the
# reference covariates.
#
# Effect estimates are reported per 1-SD increase because omics variables
# have already been standardised to z-scores.

run_single_biomarker_logistic <- function(
    data,
    outcome,
    biomarker,
    covariates
) {
  predictors <- c(covariates, biomarker)

  fit <- fit_logistic_model(
    data = data,
    outcome = outcome,
    predictors = predictors
  )

  broom::tidy(fit) %>%
    filter(term == biomarker) %>%
    mutate(
      biomarker = biomarker,
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std.error),
      CI_high = exp(estimate + 1.96 * std.error),
      p_value = p.value
    ) %>%
    select(
      biomarker,
      estimate,
      std.error,
      OR,
      CI_low,
      CI_high,
      p_value
    )
}

prevalent_biomarker_assoc_ukb <- purrr::map_dfr(
  selected_omics_main,
  ~ run_single_biomarker_logistic(
    data = ukb_bl_imp,
    outcome = "frailty_prevalent",
    biomarker = .x,
    covariates = reference_covariates
  )
) %>%
  mutate(
    cohort = "UKB",
    p_fdr = p.adjust(p_value, method = analysis_settings$p_adjust_method),
    biomarker_display = make_display_name(biomarker)
  ) %>%
  select(
    cohort,
    biomarker,
    biomarker_display,
    estimate,
    std.error,
    OR,
    CI_low,
    CI_high,
    p_value,
    p_fdr
  )

readr::write_csv(
  prevalent_biomarker_assoc_ukb,
  file.path(dir_tables, "supplementary_table3_prevalent_associations_ukb.csv")
)

# ---- 07.8 Predicted probabilities in UKB ----
# These are derivation-set predictions. External validation predictions
# are generated in Section 09 by applying UKB coefficients to ESTHER.

predict_glm_probability <- function(model, newdata) {
  stats::predict(
    model,
    newdata = newdata,
    type = "response"
  )
}

ukb_bl_imp <- ukb_bl_imp %>%
  mutate(
    pred_fi_reference = predict_glm_probability(
      trim_fi_models_ukb$Reference,
      newdata = ukb_bl_imp
    ),
    pred_fi_proteomics = predict_glm_probability(
      trim_fi_models_ukb$Reference_plus_Proteomics,
      newdata = ukb_bl_imp
    ),
    pred_fi_metabolomics = predict_glm_probability(
      trim_fi_models_ukb$Reference_plus_Metabolomics,
      newdata = ukb_bl_imp
    ),
    pred_trim_fi = predict_glm_probability(
      trim_fi_models_ukb$TRIM_FI,
      newdata = ukb_bl_imp
    )
  )

readr::write_csv(
  ukb_bl_imp,
  file.path(dir_data_processed, "ukb_baseline_imputed_z_with_trim_fi_predictions.csv")
)

# ---- 07.9 TRIM-FI coefficient table for later external validation ----
# These coefficients will be used in Section 09 to compute ESTHER
# validation predictions without refitting.

trim_fi_beta_for_validation <- broom::tidy(trim_fi_models_ukb$TRIM_FI) %>%
  select(term, estimate) %>%
  mutate(
    model = "TRIM-FI",
    .before = 1
  )

readr::write_csv(
  trim_fi_beta_for_validation,
  file.path(dir_models, "trim_fi_ukb_coefficients_for_validation.csv")
)



# ============================================================
# 08. TRIM-FP: incident frailty prediction
# ============================================================

# ---- 08.1 Analysis objective ----
# Develop the TRIM-FP model for incident frailty prediction in UKB.
#
# Longitudinal modelling:
#   - Baseline frail participants are excluded.
#   - Incident frailty is defined as first onset of frailty among
#     participants who were non-frail or pre-frail at baseline.
#   - Frailty status was assessed at discrete visits.
#   - Exact onset time was not observed.
#   - Time to incident frailty is therefore treated as interval-censored.
#   - Death is treated as a competing event.
#   - Model: interval-censored competing-risks regression.
#
# R package:
#   - intccr
#
# Function:
#   - intccr::ciregic()
#
# Expected event coding for intccr:
#   - event = 0: censored / no incident frailty
#   - event = 1: incident frailty
#   - event = 2: competing event, death
#
# Expected time interval:
#   - v: left endpoint of interval
#   - u: right endpoint of interval
#
# For incident frailty:
#   - v = last visit time at which participant was non-frail
#   - u = first visit time at which participant was frail
#
# For no incident frailty:
#   - v = last observed follow-up time
#   - u = Inf or last observed time depending on data structure accepted
#       by the analysis function.
#
# For death before frailty:
#   - event = 2
#   - interval should reflect time to death or death interval according
#     to available data.

# ---- 08.2 Follow-up data requirements ----
# The follow-up datasets should contain one row per participant or be
# transformable into one row per participant with:
#
#   id
#   follow-up FI at each visit or already-derived frailty status
#   follow-up visit time in years
#   death indicator
#   death time, if available
#
# For this manuscript:
#   - UKB frailty reassessed once at approximately 8 years post-baseline.
#   - ESTHER frailty assessed at approximately 2, 5, and 8 years.
#
# To keep the main script robust, we define a general function that accepts
# long-format follow-up data:
#
#   id
#   visit_time
#   FI
#   death
#   death_time
#
# where visit_time is years since baseline.

# ---- 08.3 Helper: classify follow-up frailty ----

classify_followup_frailty <- function(data, cohort_name) {
  data %>%
    mutate(
      frailty_fu = classify_frailty(FI, cohort_name)
    )
}

# ---- 08.4 Helper: create interval-censored competing-risk dataset ----
# Input:
#   baseline_data:
#     one row per participant, including baseline frailty status and predictors
#
#   followup_long:
#     long format, one row per participant per follow-up visit:
#       id, visit_time, FI
#
#   death_data:
#     optional columns in followup_long or separate data:
#       death, death_time
#
# Output:
#   one row per participant:
#     id, v, u, event
#
# Notes:
#   - This function assumes baseline time = 0.
#   - Participants frail at baseline are excluded.
#   - If first frail visit occurs, event = 1 and interval is
#     [last non-frail visit time, first frail visit time].
#   - If death occurs before incident frailty, event = 2.
#   - If no frailty and no death, event = 0.

make_interval_cr_data <- function(
    baseline_data,
    followup_long,
    cohort_name,
    id_var = "id",
    baseline_frailty_var = "frailty_prevalent",
    visit_time_var = "visit_time",
    FI_var = "FI",
    death_var = "death",
    death_time_var = "death_time"
) {
  # Baseline eligible participants: non-frail at baseline.
  baseline_eligible <- baseline_data %>%
    filter(.data[[baseline_frailty_var]] == 0) %>%
    select(all_of(id_var)) %>%
    distinct()

  fu <- followup_long %>%
    semi_join(baseline_eligible, by = id_var) %>%
    rename(
      visit_time_tmp = all_of(visit_time_var),
      FI_tmp = all_of(FI_var)
    ) %>%
    filter(!is.na(visit_time_tmp)) %>%
    mutate(
      frailty_fu_tmp = classify_frailty(FI_tmp, cohort_name)
    ) %>%
    arrange(.data[[id_var]], visit_time_tmp)

  # Extract death information if available.
  has_death <- all(c(death_var, death_time_var) %in% names(followup_long))

  if (has_death) {
    death_info <- followup_long %>%
      semi_join(baseline_eligible, by = id_var) %>%
      group_by(.data[[id_var]]) %>%
      summarise(
        death_tmp = as.integer(any(.data[[death_var]] == 1, na.rm = TRUE)),
        death_time_tmp = {
          death_times <- .data[[death_time_var]][.data[[death_var]] == 1]
          death_times <- death_times[!is.na(death_times)]
          if (length(death_times) > 0) {
            min(death_times)
          } else {
            NA_real_
          }
        },
        .groups = "drop"
      )
  } else {
    death_info <- baseline_eligible %>%
      mutate(
        death_tmp = 0L,
        death_time_tmp = NA_real_
      )
  }

  fu_summary <- fu %>%
    group_by(.data[[id_var]]) %>%
    summarise(
      any_frail = any(frailty_fu_tmp == 1, na.rm = TRUE),
      first_frail_time = {
        frail_times <- visit_time_tmp[frailty_fu_tmp == 1]
        frail_times <- frail_times[!is.na(frail_times)]
        if (length(frail_times) > 0) {
          min(frail_times)
        } else {
          NA_real_
        }
      },
      last_observed_time = {
        observed_times <- visit_time_tmp[!is.na(visit_time_tmp)]
        if (length(observed_times) > 0) {
          max(observed_times)
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    )

  last_nonfrail <- fu %>%
    left_join(
      fu_summary %>%
        select(all_of(id_var), first_frail_time),
      by = id_var
    ) %>%
    filter(
      !is.na(first_frail_time),
      frailty_fu_tmp == 0,
      visit_time_tmp < first_frail_time
    ) %>%
    group_by(.data[[id_var]]) %>%
    summarise(
      last_nonfrail_before_frail = max(visit_time_tmp, na.rm = TRUE),
      .groups = "drop"
    )

  interval_data <- fu_summary %>%
    left_join(last_nonfrail, by = id_var) %>%
    mutate(
      last_nonfrail_before_frail = ifelse(
        is.na(last_nonfrail_before_frail),
        0,
        last_nonfrail_before_frail
      )
    ) %>%
    left_join(death_info, by = id_var) %>%
    mutate(
      death_tmp = ifelse(is.na(death_tmp), 0L, death_tmp),

      event = case_when(
        any_frail & (is.na(death_time_tmp) | first_frail_time <= death_time_tmp) ~ 1L,
        any_frail & death_tmp == 1L & death_time_tmp < first_frail_time ~ 2L,
        !any_frail & death_tmp == 1L ~ 2L,
        TRUE ~ 0L
      ),

      v = case_when(
        event == 1L ~ last_nonfrail_before_frail,
        event == 2L ~ death_time_tmp,
        event == 0L ~ last_observed_time,
        TRUE ~ NA_real_
      ),

      u = case_when(
        event == 1L ~ first_frail_time,
        event == 2L ~ death_time_tmp,
        event == 0L ~ last_observed_time,
        TRUE ~ NA_real_
      )
    ) %>%
    select(
      all_of(id_var),
      v,
      u,
      event,
      first_frail_time,
      last_observed_time,
      death_tmp,
      death_time_tmp
    )

  if (any(interval_data$event == 1 & interval_data$v > interval_data$u, na.rm = TRUE)) {
    stop("Invalid incident frailty intervals detected: v > u.")
  }

  if (any(interval_data$event == 1 & interval_data$v == interval_data$u, na.rm = TRUE)) {
    warning("Some incident frailty intervals have v == u. Check follow-up timing.")
  }

  interval_data
}

# ---- 08.5 Prepare follow-up long datasets ----
# This section assumes ukb_fu and esther_fu are already in long format with:
#   id, visit_time, FI, death, death_time
#
# If your current follow-up data are wide format, reshape them before this step.
#
# Example wide-to-long template:
#
# esther_fu_long <- esther_fu %>%
#   pivot_longer(
#     cols = starts_with("FI_y"),
#     names_to = "visit",
#     values_to = "FI"
#   ) %>%
#   mutate(
#     visit_time = case_when(
#       visit == "FI_y2" ~ 2,
#       visit == "FI_y5" ~ 5,
#       visit == "FI_y8" ~ 8,
#       TRUE ~ NA_real_
#     )
#   )
#
# ukb_fu_long <- ukb_fu %>%
#   transmute(
#     id = id,
#     visit_time = 8,
#     FI = FI_fu,
#     death = death,
#     death_time = death_time
#   )

ukb_fu_long <- ukb_fu
esther_fu_long <- esther_fu

# ---- 08.6 Construct interval-censored competing-risk datasets ----

ukb_interval <- make_interval_cr_data(
  baseline_data = ukb_bl_imp,
  followup_long = ukb_fu_long,
  cohort_name = "UKB",
  id_var = "id",
  baseline_frailty_var = "frailty_prevalent",
  visit_time_var = "visit_time",
  FI_var = "FI",
  death_var = "death",
  death_time_var = "death_time"
)

esther_interval <- make_interval_cr_data(
  baseline_data = esther_bl_imp,
  followup_long = esther_fu_long,
  cohort_name = "ESTHER",
  id_var = "id",
  baseline_frailty_var = "frailty_prevalent",
  visit_time_var = "visit_time",
  FI_var = "FI",
  death_var = "death",
  death_time_var = "death_time"
)

# ---- 08.7 Merge baseline predictors with interval outcomes ----

ukb_fp_data <- ukb_bl_imp %>%
  inner_join(ukb_interval, by = "id") %>%
  filter(frailty_prevalent == 0)

esther_fp_data <- esther_bl_imp %>%
  inner_join(esther_interval, by = "id") %>%
  filter(frailty_prevalent == 0)

# ---- 08.8 Follow-up analysis check ----
# Expected manuscript numbers after inclusion/exclusion:
#   UKB:    n = 2,546; incident frailty = 250; 9.8%
#   ESTHER: n = 2,593; incident frailty = 255; 9.8%

incident_check <- bind_rows(
  ukb_fp_data %>%
    summarise(
      cohort = "UKB",
      n = n(),
      incident_frailty_n = sum(event == 1, na.rm = TRUE),
      death_competing_n = sum(event == 2, na.rm = TRUE),
      censored_n = sum(event == 0, na.rm = TRUE),
      incident_percent = 100 * mean(event == 1, na.rm = TRUE),
      median_followup = median(u, na.rm = TRUE)
    ),

  esther_fp_data %>%
    summarise(
      cohort = "ESTHER",
      n = n(),
      incident_frailty_n = sum(event == 1, na.rm = TRUE),
      death_competing_n = sum(event == 2, na.rm = TRUE),
      censored_n = sum(event == 0, na.rm = TRUE),
      incident_percent = 100 * mean(event == 1, na.rm = TRUE),
      median_followup = median(u, na.rm = TRUE)
    )
)

print(incident_check)

readr::write_csv(
  incident_check,
  file.path(dir_logs, "incident_frailty_check.csv")
)

# ---- 08.9 Define TRIM-FP model variable sets ----

trim_fp_variable_sets <- list(
  Reference = reference_covariates,
  Reference_plus_Proteomics = c(reference_covariates, main_analysis$proteins),
  Reference_plus_Metabolomics = c(reference_covariates, main_analysis$metabolites),
  TRIM_FP = c(reference_covariates, main_analysis$omics)
)

# ---- 08.10 Helper: fit interval-censored competing-risk model ----
# intccr::ciregic() uses Surv2(v, u, event) for interval-censored
# competing-risks data.
#
# alpha = c(0, 0) specifies the model form used in the package examples
# for semiparametric proportional subdistribution hazards regression.

fit_ciregic_model <- function(data, predictors) {
  vars_needed <- unique(c("v", "u", "event", predictors))

  dat <- data %>%
    select(all_of(vars_needed)) %>%
    filter(
      !is.na(v),
      !is.na(u),
      !is.na(event),
      u > v
    )

  dat <- apply_category_levels(dat, category_levels)

  formula <- as.formula(
    paste(
      "intccr::Surv2(v, u, event) ~",
      paste(predictors, collapse = " + ")
    )
  )

  intccr::ciregic(
    formula = formula,
    data = dat,
    alpha = c(0, 0)
  )
}

# ---- 08.11 Fit TRIM-FP models in UKB ----

trim_fp_models_ukb <- list(
  Reference = fit_ciregic_model(
    data = ukb_fp_data,
    predictors = trim_fp_variable_sets$Reference
  ),

  Reference_plus_Proteomics = fit_ciregic_model(
    data = ukb_fp_data,
    predictors = trim_fp_variable_sets$Reference_plus_Proteomics
  ),

  Reference_plus_Metabolomics = fit_ciregic_model(
    data = ukb_fp_data,
    predictors = trim_fp_variable_sets$Reference_plus_Metabolomics
  ),

  TRIM_FP = fit_ciregic_model(
    data = ukb_fp_data,
    predictors = trim_fp_variable_sets$TRIM_FP
  )
)

saveRDS(
  trim_fp_models_ukb,
  file.path(dir_models, "trim_fp_models_ukb.rds")
)

# ---- 08.12 Extract TRIM-FP coefficients ----
# ciregic object structure may vary by package version.
# This helper attempts to extract coefficients robustly.

extract_ciregic_coefficients <- function(model, model_name) {
  beta <- stats::coef(model)

  tibble(
    model = model_name,
    term = names(beta),
    estimate = as.numeric(beta),
    sHR = exp(estimate)
  )
}

trim_fp_coefficients_ukb <- purrr::imap_dfr(
  trim_fp_models_ukb,
  extract_ciregic_coefficients
)

readr::write_csv(
  trim_fp_coefficients_ukb,
  file.path(dir_tables, "trim_fp_coefficients_ukb.csv")
)

# ---- 08.13 Single-biomarker incident frailty associations ----
# For Figure 2 / Supplementary Table 4:
# each selected biomarker is modelled separately, adjusted for reference
# covariates, using interval-censored competing-risks regression.

run_single_biomarker_ciregic <- function(
    data,
    biomarker,
    covariates
) {
  predictors <- c(covariates, biomarker)

  fit <- fit_ciregic_model(
    data = data,
    predictors = predictors
  )

  beta <- stats::coef(fit)

  if (!biomarker %in% names(beta)) {
    return(
      tibble(
        biomarker = biomarker,
        estimate = NA_real_,
        sHR = NA_real_
      )
    )
  }

  tibble(
    biomarker = biomarker,
    estimate = as.numeric(beta[biomarker]),
    sHR = exp(estimate)
  )
}

incident_biomarker_assoc_ukb <- purrr::map_dfr(
  selected_omics_main,
  ~ run_single_biomarker_ciregic(
    data = ukb_fp_data,
    biomarker = .x,
    covariates = reference_covariates
  )
) %>%
  mutate(
    cohort = "UKB",
    biomarker_display = make_display_name(biomarker)
  ) %>%
  select(
    cohort,
    biomarker,
    biomarker_display,
    estimate,
    sHR
  )

readr::write_csv(
  incident_biomarker_assoc_ukb,
  file.path(dir_tables, "supplementary_table4_incident_associations_ukb.csv")
)

# ---- 08.14 Linear predictors / relative risk scores in UKB ----
# For TRIM-FP, the model provides a relative risk score.
# Absolute 8-year risk requires baseline CIF / calibration.

make_model_matrix_for_predictors <- function(data, predictors) {
  dat <- data %>%
    select(all_of(predictors))

  dat <- apply_category_levels(dat, category_levels)

  formula_x <- as.formula(
    paste("~", paste(predictors, collapse = " + "))
  )

  model.matrix(formula_x, data = dat)[, -1, drop = FALSE]
}

linear_predictor_from_beta <- function(data, predictors, beta) {
  x <- make_model_matrix_for_predictors(data, predictors)

  common_terms <- intersect(colnames(x), names(beta))

  lp <- as.numeric(x[, common_terms, drop = FALSE] %*% beta[common_terms])

  # Add intercept if present.
  if ("(Intercept)" %in% names(beta)) {
    lp <- lp + beta["(Intercept)"]
  }

  lp
}

# For ciregic models, use coefficient names directly.
beta_fp_reference <- stats::coef(trim_fp_models_ukb$Reference)
beta_fp_proteomics <- stats::coef(trim_fp_models_ukb$Reference_plus_Proteomics)
beta_fp_metabolomics <- stats::coef(trim_fp_models_ukb$Reference_plus_Metabolomics)
beta_fp_trim <- stats::coef(trim_fp_models_ukb$TRIM_FP)

ukb_fp_data <- ukb_fp_data %>%
  mutate(
    lp_fp_reference = linear_predictor_from_beta(
      data = ukb_fp_data,
      predictors = trim_fp_variable_sets$Reference,
      beta = beta_fp_reference
    ),
    lp_fp_proteomics = linear_predictor_from_beta(
      data = ukb_fp_data,
      predictors = trim_fp_variable_sets$Reference_plus_Proteomics,
      beta = beta_fp_proteomics
    ),
    lp_fp_metabolomics = linear_predictor_from_beta(
      data = ukb_fp_data,
      predictors = trim_fp_variable_sets$Reference_plus_Metabolomics,
      beta = beta_fp_metabolomics
    ),
    lp_trim_fp = linear_predictor_from_beta(
      data = ukb_fp_data,
      predictors = trim_fp_variable_sets$TRIM_FP,
      beta = beta_fp_trim
    ),
    rr_trim_fp = exp(lp_trim_fp)
  )

readr::write_csv(
  ukb_fp_data,
  file.path(dir_data_processed, "ukb_followup_interval_with_trim_fp_scores.csv")
)

# ---- 08.15 TRIM-FP coefficient table for external validation ----

trim_fp_beta_for_validation <- tibble(
  model = "TRIM-FP",
  term = names(beta_fp_trim),
  estimate = as.numeric(beta_fp_trim)
)

readr::write_csv(
  trim_fp_beta_for_validation,
  file.path(dir_models, "trim_fp_ukb_coefficients_for_validation.csv")
)


# ============================================================
# 09. External validation in ESTHER
# ============================================================

# ---- 09.1 Validation principle ----
# External validation:
#
#   - Regression coefficients estimated in UKB are directly applied
#     to ESTHER.
#   - The main TRIM models are not refitted in ESTHER.
#
# This applies to:
#   - TRIM-FI for prevalent frailty identification
#   - TRIM-FP for incident frailty prediction

# ---- 09.2 Ensure ESTHER categorical coding matches UKB ----

esther_bl_imp <- apply_category_levels(esther_bl_imp, category_levels)
esther_fp_data <- apply_category_levels(esther_fp_data, category_levels)

# ---- 09.3 Apply UKB TRIM-FI coefficients to ESTHER ----
# For logistic regression, predict.glm with newdata uses the UKB-fitted
# model coefficients directly and returns ESTHER validation probabilities.

esther_bl_imp <- esther_bl_imp %>%
  mutate(
    pred_fi_reference = predict_glm_probability(
      trim_fi_models_ukb$Reference,
      newdata = esther_bl_imp
    ),
    pred_fi_proteomics = predict_glm_probability(
      trim_fi_models_ukb$Reference_plus_Proteomics,
      newdata = esther_bl_imp
    ),
    pred_fi_metabolomics = predict_glm_probability(
      trim_fi_models_ukb$Reference_plus_Metabolomics,
      newdata = esther_bl_imp
    ),
    pred_trim_fi = predict_glm_probability(
      trim_fi_models_ukb$TRIM_FI,
      newdata = esther_bl_imp
    )
  )

readr::write_csv(
  esther_bl_imp,
  file.path(dir_data_processed, "esther_baseline_imputed_z_with_trim_fi_predictions.csv")
)

# ---- 09.4 Apply UKB TRIM-FP coefficients to ESTHER ----
# For TRIM-FP, we apply UKB ciregic coefficients to ESTHER predictors
# to obtain linear predictors and relative risk scores.
#
# Absolute risk prediction requires baseline CIF / calibration and is
# handled in Section 10 if baseline risk information is available.

esther_fp_data <- esther_fp_data %>%
  mutate(
    lp_fp_reference = linear_predictor_from_beta(
      data = esther_fp_data,
      predictors = trim_fp_variable_sets$Reference,
      beta = beta_fp_reference
    ),
    lp_fp_proteomics = linear_predictor_from_beta(
      data = esther_fp_data,
      predictors = trim_fp_variable_sets$Reference_plus_Proteomics,
      beta = beta_fp_proteomics
    ),
    lp_fp_metabolomics = linear_predictor_from_beta(
      data = esther_fp_data,
      predictors = trim_fp_variable_sets$Reference_plus_Metabolomics,
      beta = beta_fp_metabolomics
    ),
    lp_trim_fp = linear_predictor_from_beta(
      data = esther_fp_data,
      predictors = trim_fp_variable_sets$TRIM_FP,
      beta = beta_fp_trim
    ),
    rr_trim_fp = exp(lp_trim_fp)
  )

readr::write_csv(
  esther_fp_data,
  file.path(dir_data_processed, "esther_followup_interval_with_trim_fp_scores.csv")
)

# ---- 09.5 Single-biomarker associations in ESTHER ----
# These are not the external validation of the full TRIM model.
# They are cohort-specific association estimates used for Figure 2 and
# Supplementary Tables 3 and 4.

prevalent_biomarker_assoc_esther <- purrr::map_dfr(
  selected_omics_main,
  ~ run_single_biomarker_logistic(
    data = esther_bl_imp,
    outcome = "frailty_prevalent",
    biomarker = .x,
    covariates = reference_covariates
  )
) %>%
  mutate(
    cohort = "ESTHER",
    p_fdr = p.adjust(p_value, method = analysis_settings$p_adjust_method),
    biomarker_display = make_display_name(biomarker)
  ) %>%
  select(
    cohort,
    biomarker,
    biomarker_display,
    estimate,
    std.error,
    OR,
    CI_low,
    CI_high,
    p_value,
    p_fdr
  )

readr::write_csv(
  prevalent_biomarker_assoc_esther,
  file.path(dir_tables, "supplementary_table3_prevalent_associations_esther.csv")
)

incident_biomarker_assoc_esther <- purrr::map_dfr(
  selected_omics_main,
  ~ run_single_biomarker_ciregic(
    data = esther_fp_data,
    biomarker = .x,
    covariates = reference_covariates
  )
) %>%
  mutate(
    cohort = "ESTHER",
    biomarker_display = make_display_name(biomarker)
  ) %>%
  select(
    cohort,
    biomarker,
    biomarker_display,
    estimate,
    sHR
  )

readr::write_csv(
  incident_biomarker_assoc_esther,
  file.path(dir_tables, "supplementary_table4_incident_associations_esther.csv")
)

# ---- 09.6 Combine biomarker association tables ----

prevalent_biomarker_assoc_all <- bind_rows(
  prevalent_biomarker_assoc_ukb,
  prevalent_biomarker_assoc_esther
)

incident_biomarker_assoc_all <- bind_rows(
  incident_biomarker_assoc_ukb,
  incident_biomarker_assoc_esther
)

readr::write_csv(
  prevalent_biomarker_assoc_all,
  file.path(dir_tables, "supplementary_table3_prevalent_associations_all.csv")
)

readr::write_csv(
  incident_biomarker_assoc_all,
  file.path(dir_tables, "supplementary_table4_incident_associations_all.csv")
)

# ---- 09.7 Validation prediction checks ----

validation_prediction_check <- bind_rows(
  esther_bl_imp %>%
    summarise(
      cohort = "ESTHER",
      outcome = "Prevalent frailty",
      n = n(),
      event_n = sum(frailty_prevalent == 1, na.rm = TRUE),
      pred_reference_mean = mean(pred_fi_reference, na.rm = TRUE),
      pred_trim_mean = mean(pred_trim_fi, na.rm = TRUE),
      pred_reference_sd = sd(pred_fi_reference, na.rm = TRUE),
      pred_trim_sd = sd(pred_trim_fi, na.rm = TRUE)
    ),

  esther_fp_data %>%
    summarise(
      cohort = "ESTHER",
      outcome = "Incident frailty",
      n = n(),
      event_n = sum(event == 1, na.rm = TRUE),
      pred_reference_mean = mean(lp_fp_reference, na.rm = TRUE),
      pred_trim_mean = mean(lp_trim_fp, na.rm = TRUE),
      pred_reference_sd = sd(lp_fp_reference, na.rm = TRUE),
      pred_trim_sd = sd(lp_trim_fp, na.rm = TRUE)
    )
)

print(validation_prediction_check)

readr::write_csv(
  validation_prediction_check,
  file.path(dir_logs, "esther_external_validation_prediction_check.csv")
)


# ============================================================
# 10. Model performance
# ============================================================

# ---- 10.1 Analysis objective ----
# Compare model performance between:
#   1. Reference model
#   2. Reference + Proteomics
#   3. Reference + Metabolomics
#   4. TRIM model
#
# For prevalent frailty:
#   - AUC
#   - DeLong test for AUC difference
#   - NRI / IDI
#   - calibration
#   - decision curve analysis
#
# For incident frailty:
#   - C-statistic / time-dependent discrimination
#   - bootstrap test for C-statistic difference
#
# Incident frailty NRI, IDI, calibration, and DCA require calibrated
# absolute 8-year risks and are not computed from uncalibrated linear
# predictors in this public script.
#
# Important:
#   - Baseline frailty models use predicted probabilities.
#   - Incident frailty models use linear predictors / risk scores unless
#     calibrated 8-year absolute risks are available.

# ---- 10.2 Helper: AUC for prevalent frailty ----

calc_auc_binary <- function(data, outcome, prediction, cohort, model) {
  roc_obj <- pROC::roc(
    response = data[[outcome]],
    predictor = data[[prediction]],
    quiet = TRUE,
    direction = "<"
  )

  ci_obj <- pROC::ci.auc(roc_obj)

  tibble(
    cohort = cohort,
    outcome = outcome,
    model = model,
    metric = "AUC",
    estimate = as.numeric(pROC::auc(roc_obj)),
    CI_low = as.numeric(ci_obj[1]),
    CI_high = as.numeric(ci_obj[3])
  )
}

# ---- 10.3 AUC for TRIM-FI in UKB and ESTHER ----

auc_fi_results <- bind_rows(
  calc_auc_binary(ukb_bl_imp, "frailty_prevalent", "pred_fi_reference", "UKB", "Reference"),
  calc_auc_binary(ukb_bl_imp, "frailty_prevalent", "pred_fi_proteomics", "UKB", "Reference + Proteomics"),
  calc_auc_binary(ukb_bl_imp, "frailty_prevalent", "pred_fi_metabolomics", "UKB", "Reference + Metabolomics"),
  calc_auc_binary(ukb_bl_imp, "frailty_prevalent", "pred_trim_fi", "UKB", "TRIM-FI"),

  calc_auc_binary(esther_bl_imp, "frailty_prevalent", "pred_fi_reference", "ESTHER", "Reference"),
  calc_auc_binary(esther_bl_imp, "frailty_prevalent", "pred_fi_proteomics", "ESTHER", "Reference + Proteomics"),
  calc_auc_binary(esther_bl_imp, "frailty_prevalent", "pred_fi_metabolomics", "ESTHER", "Reference + Metabolomics"),
  calc_auc_binary(esther_bl_imp, "frailty_prevalent", "pred_trim_fi", "ESTHER", "TRIM-FI")
)

readr::write_csv(
  auc_fi_results,
  file.path(dir_tables, "table2_auc_prevalent_frailty.csv")
)

# ---- 10.4 DeLong tests for AUC differences ----

compare_auc_delong <- function(data, outcome, pred_reference, pred_new, cohort, model_new) {
  roc_ref <- pROC::roc(
    response = data[[outcome]],
    predictor = data[[pred_reference]],
    quiet = TRUE,
    direction = "<"
  )

  roc_new <- pROC::roc(
    response = data[[outcome]],
    predictor = data[[pred_new]],
    quiet = TRUE,
    direction = "<"
  )

  test_obj <- pROC::roc.test(
    roc_ref,
    roc_new,
    method = "delong",
    paired = TRUE
  )

  tibble(
    cohort = cohort,
    outcome = outcome,
    comparison = paste(model_new, "vs Reference"),
    auc_reference = as.numeric(pROC::auc(roc_ref)),
    auc_new = as.numeric(pROC::auc(roc_new)),
    delta_auc = auc_new - auc_reference,
    p_value = as.numeric(test_obj$p.value)
  )
}

auc_fi_differences <- bind_rows(
  compare_auc_delong(ukb_bl_imp, "frailty_prevalent", "pred_fi_reference", "pred_fi_proteomics", "UKB", "Reference + Proteomics"),
  compare_auc_delong(ukb_bl_imp, "frailty_prevalent", "pred_fi_reference", "pred_fi_metabolomics", "UKB", "Reference + Metabolomics"),
  compare_auc_delong(ukb_bl_imp, "frailty_prevalent", "pred_fi_reference", "pred_trim_fi", "UKB", "TRIM-FI"),

  compare_auc_delong(esther_bl_imp, "frailty_prevalent", "pred_fi_reference", "pred_fi_proteomics", "ESTHER", "Reference + Proteomics"),
  compare_auc_delong(esther_bl_imp, "frailty_prevalent", "pred_fi_reference", "pred_fi_metabolomics", "ESTHER", "Reference + Metabolomics"),
  compare_auc_delong(esther_bl_imp, "frailty_prevalent", "pred_fi_reference", "pred_trim_fi", "ESTHER", "TRIM-FI")
)

readr::write_csv(
  auc_fi_differences,
  file.path(dir_tables, "table2_delta_auc_prevalent_frailty.csv")
)

# ---- 10.5 Incident frailty discrimination ----
# The manuscript reports time-dependent Harrell's C-statistic for
# incident frailty.
#
# Public implementation note:
# This helper provides an approximate discrimination estimate using
# survival::concordance on a time-to-first-event approximation:
#   - time = u
#   - status = 1 for incident frailty
#   - death and censoring are treated as non-frailty events for
#     discrimination of incident frailty.
#
# For replication, replace this helper with the
# original time-dependent C-statistic implementation based on timeROC.

calc_cstat_incident <- function(data, score, cohort, model) {
  dat <- data %>%
    filter(!is.na(u), !is.na(event), !is.na(.data[[score]])) %>%
    mutate(
      incident_status = as.integer(event == 1),
      risk_score = .data[[score]]
    )

  concordance_obj <- survival::concordance(
    survival::Surv(u, incident_status) ~ risk_score,
    data = dat,
    reverse = TRUE
  )

  tibble(
    cohort = cohort,
    outcome = "incident_frailty",
    model = model,
    metric = "C-statistic",
    estimate = as.numeric(concordance_obj$concordance)
  )
}

cstat_fp_results <- bind_rows(
  calc_cstat_incident(ukb_fp_data, "lp_fp_reference", "UKB", "Reference"),
  calc_cstat_incident(ukb_fp_data, "lp_fp_proteomics", "UKB", "Reference + Proteomics"),
  calc_cstat_incident(ukb_fp_data, "lp_fp_metabolomics", "UKB", "Reference + Metabolomics"),
  calc_cstat_incident(ukb_fp_data, "lp_trim_fp", "UKB", "TRIM-FP"),

  calc_cstat_incident(esther_fp_data, "lp_fp_reference", "ESTHER", "Reference"),
  calc_cstat_incident(esther_fp_data, "lp_fp_proteomics", "ESTHER", "Reference + Proteomics"),
  calc_cstat_incident(esther_fp_data, "lp_fp_metabolomics", "ESTHER", "Reference + Metabolomics"),
  calc_cstat_incident(esther_fp_data, "lp_trim_fp", "ESTHER", "TRIM-FP")
)

readr::write_csv(
  cstat_fp_results,
  file.path(dir_tables, "table2_cstat_incident_frailty.csv")
)

# ---- 10.6 Bootstrap difference in incident C-statistics ----

bootstrap_delta_cstat <- function(
    data,
    score_reference,
    score_new,
    B = 1000,
    seed = 20260430
) {
  set.seed(seed)

  n <- nrow(data)
  delta <- numeric(B)

  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    dat_b <- data[idx, ]

    c_ref <- calc_cstat_incident(dat_b, score_reference, "bootstrap", "Reference")$estimate
    c_new <- calc_cstat_incident(dat_b, score_new, "bootstrap", "New")$estimate

    delta[b] <- c_new - c_ref
  }

  tibble(
    delta = mean(delta, na.rm = TRUE),
    CI_low = quantile(delta, 0.025, na.rm = TRUE),
    CI_high = quantile(delta, 0.975, na.rm = TRUE),
    p_value = 2 * min(
      mean(delta <= 0, na.rm = TRUE),
      mean(delta >= 0, na.rm = TRUE)
    )
  )
}

cstat_fp_differences <- bind_rows(
  bootstrap_delta_cstat(ukb_fp_data, "lp_fp_reference", "lp_fp_proteomics") %>%
    mutate(cohort = "UKB", comparison = "Reference + Proteomics vs Reference"),
  bootstrap_delta_cstat(ukb_fp_data, "lp_fp_reference", "lp_fp_metabolomics") %>%
    mutate(cohort = "UKB", comparison = "Reference + Metabolomics vs Reference"),
  bootstrap_delta_cstat(ukb_fp_data, "lp_fp_reference", "lp_trim_fp") %>%
    mutate(cohort = "UKB", comparison = "TRIM-FP vs Reference"),

  bootstrap_delta_cstat(esther_fp_data, "lp_fp_reference", "lp_fp_proteomics") %>%
    mutate(cohort = "ESTHER", comparison = "Reference + Proteomics vs Reference"),
  bootstrap_delta_cstat(esther_fp_data, "lp_fp_reference", "lp_fp_metabolomics") %>%
    mutate(cohort = "ESTHER", comparison = "Reference + Metabolomics vs Reference"),
  bootstrap_delta_cstat(esther_fp_data, "lp_fp_reference", "lp_trim_fp") %>%
    mutate(cohort = "ESTHER", comparison = "TRIM-FP vs Reference")
)

readr::write_csv(
  cstat_fp_differences,
  file.path(dir_tables, "table2_delta_cstat_incident_frailty.csv")
)

# ---- 10.7 Categorical NRI helper ----
# Risk categories:
#   <5%, 5-10%, >=10%
#
# This helper calculates category-based net reclassification improvement
# for binary outcomes using absolute predicted risks.
#
# Intended use:
#   - TRIM-FI: predicted probabilities from logistic regression.
#   - TRIM-FP: calibrated absolute risks only, if available.
#
# Not intended use:
#   - Do not apply this function to uncalibrated linear predictors,
#     relative risk scores, log hazards, or arbitrary risk scores.
#
# Definitions:
#   For participants with events:
#     NRI_events = P(up-classified | event) - P(down-classified | event)
#
#   For participants without events:
#     NRI_nonevents = P(down-classified | nonevent) -
#                     P(up-classified | nonevent)
#
#   Overall categorical NRI:
#     NRI_total = NRI_events + NRI_nonevents

risk_category <- function(
    p,
    cutoffs = analysis_settings$nri_cutoffs,
    labels = NULL
) {
  if (!is.numeric(p)) {
    stop("`p` must be a numeric vector of absolute predicted risks.")
  }

  if (!is.numeric(cutoffs) || length(cutoffs) < 1) {
    stop("`cutoffs` must be a numeric vector with at least one cutoff.")
  }

  if (any(is.na(cutoffs))) {
    stop("`cutoffs` must not contain missing values.")
  }

  if (any(cutoffs <= 0 | cutoffs >= 1)) {
    stop("All `cutoffs` must be between 0 and 1.")
  }

  if (is.unsorted(cutoffs, strictly = TRUE)) {
    stop("`cutoffs` must be strictly increasing.")
  }

  if (any(!is.na(p) & (p < 0 | p > 1))) {
    stop(
      "`p` contains values outside [0, 1]. ",
      "Categorical NRI requires absolute predicted risks or probabilities."
    )
  }

  if (is.null(labels)) {
    labels <- c(
      paste0("<", scales::percent(cutoffs[1], accuracy = 1)),
      paste0(
        scales::percent(cutoffs[-length(cutoffs)], accuracy = 1),
        "-<",
        scales::percent(cutoffs[-1], accuracy = 1)
      ),
      paste0(">=", scales::percent(cutoffs[length(cutoffs)], accuracy = 1))
    )
  }

  expected_n_labels <- length(cutoffs) + 1

  if (length(labels) != expected_n_labels) {
    stop(
      "`labels` must have length equal to length(cutoffs) + 1."
    )
  }

  cut(
    p,
    breaks = c(-Inf, cutoffs, Inf),
    labels = labels,
    right = FALSE,
    ordered_result = TRUE
  )
}

calc_categorical_nri <- function(
    y,
    p_ref,
    p_new,
    cutoffs = analysis_settings$nri_cutoffs,
    labels = NULL
) {
  if (length(y) != length(p_ref) || length(y) != length(p_new)) {
    stop("`y`, `p_ref`, and `p_new` must have the same length.")
  }

  if (!all(na.omit(y) %in% c(0, 1))) {
    stop("`y` must be coded as 0/1.")
  }

  if (!is.numeric(p_ref) || !is.numeric(p_new)) {
    stop("`p_ref` and `p_new` must be numeric predicted risks.")
  }

  if (any(!is.na(p_ref) & (p_ref < 0 | p_ref > 1)) ||
      any(!is.na(p_new) & (p_new < 0 | p_new > 1))) {
    stop(
      "`p_ref` and `p_new` must be absolute predicted risks in [0, 1]. ",
      "Do not use linear predictors or relative risk scores."
    )
  }

  dat <- tibble(
    y = as.integer(y),
    p_ref = p_ref,
    p_new = p_new
  ) %>%
    filter(
      !is.na(y),
      !is.na(p_ref),
      !is.na(p_new)
    )

  if (nrow(dat) == 0) {
    stop("No complete observations are available for NRI calculation.")
  }

  if (!any(dat$y == 1)) {
    stop("No events are available for NRI calculation.")
  }

  if (!any(dat$y == 0)) {
    stop("No non-events are available for NRI calculation.")
  }

  dat <- dat %>%
    mutate(
      category_ref = risk_category(
        p = p_ref,
        cutoffs = cutoffs,
        labels = labels
      ),
      category_new = risk_category(
        p = p_new,
        cutoffs = cutoffs,
        labels = labels
      ),
      category_ref_num = as.integer(category_ref),
      category_new_num = as.integer(category_new),
      reclassification = case_when(
        category_new_num > category_ref_num ~ "Up",
        category_new_num < category_ref_num ~ "Down",
        category_new_num == category_ref_num ~ "Unchanged",
        TRUE ~ NA_character_
      )
    )

  events <- dat$y == 1
  nonevents <- dat$y == 0

  event_up <- mean(dat$category_new_num[events] > dat$category_ref_num[events])
  event_down <- mean(dat$category_new_num[events] < dat$category_ref_num[events])
  event_unchanged <- mean(dat$category_new_num[events] == dat$category_ref_num[events])

  nonevent_up <- mean(dat$category_new_num[nonevents] > dat$category_ref_num[nonevents])
  nonevent_down <- mean(dat$category_new_num[nonevents] < dat$category_ref_num[nonevents])
  nonevent_unchanged <- mean(dat$category_new_num[nonevents] == dat$category_ref_num[nonevents])

  nri_events <- event_up - event_down
  nri_nonevents <- nonevent_down - nonevent_up
  nri_total <- nri_events + nri_nonevents

  tibble(
    n_total = nrow(dat),
    n_events = sum(events),
    n_nonevents = sum(nonevents),

    event_up = event_up,
    event_down = event_down,
    event_unchanged = event_unchanged,

    nonevent_up = nonevent_up,
    nonevent_down = nonevent_down,
    nonevent_unchanged = nonevent_unchanged,

    nri_events = nri_events,
    nri_nonevents = nri_nonevents,
    categorical_nri = nri_total
  )
}

# ---- 10.8 IDI helper for binary outcomes ----

calc_idi_binary <- function(y, p_ref, p_new) {
  y <- as.integer(y)

  events <- y == 1
  nonevents <- y == 0

  discrimination_ref <- mean(p_ref[events], na.rm = TRUE) -
    mean(p_ref[nonevents], na.rm = TRUE)

  discrimination_new <- mean(p_new[events], na.rm = TRUE) -
    mean(p_new[nonevents], na.rm = TRUE)

  tibble(
    discrimination_ref = discrimination_ref,
    discrimination_new = discrimination_new,
    idi = discrimination_new - discrimination_ref
  )
}

# ---- 10.9 Reclassification metrics for TRIM-FI ----

reclassification_fi <- bind_rows(
  calc_categorical_nri(
    y = ukb_bl_imp$frailty_prevalent,
    p_ref = ukb_bl_imp$pred_fi_reference,
    p_new = ukb_bl_imp$pred_trim_fi
  ) %>%
    bind_cols(
      calc_idi_binary(
        y = ukb_bl_imp$frailty_prevalent,
        p_ref = ukb_bl_imp$pred_fi_reference,
        p_new = ukb_bl_imp$pred_trim_fi
      )
    ) %>%
    mutate(cohort = "UKB", outcome = "Prevalent frailty", model = "TRIM-FI"),

  calc_categorical_nri(
    y = esther_bl_imp$frailty_prevalent,
    p_ref = esther_bl_imp$pred_fi_reference,
    p_new = esther_bl_imp$pred_trim_fi
  ) %>%
    bind_cols(
      calc_idi_binary(
        y = esther_bl_imp$frailty_prevalent,
        p_ref = esther_bl_imp$pred_fi_reference,
        p_new = esther_bl_imp$pred_trim_fi
      )
    ) %>%
    mutate(cohort = "ESTHER", outcome = "Prevalent frailty", model = "TRIM-FI")
)

readr::write_csv(
  reclassification_fi,
  file.path(dir_tables, "reclassification_prevalent_frailty_trim_fi.csv")
)

# ---- 10.10 Calibration plot data ----
# For baseline frailty, observed outcome is binary and predicted risk is
# logistic probability.
#
# For incident frailty, calibration should use calibrated 8-year risk.
# If calibrated risk is not available, this function should not be used
# for TRIM-FP manuscript calibration.

make_calibration_data_binary <- function(data, outcome, prediction, n_groups = 10) {
  data %>%
    filter(!is.na(.data[[outcome]]), !is.na(.data[[prediction]])) %>%
    mutate(
      risk_group = dplyr::ntile(.data[[prediction]], n_groups)
    ) %>%
    group_by(risk_group) %>%
    summarise(
      predicted = mean(.data[[prediction]], na.rm = TRUE),
      observed = mean(.data[[outcome]], na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
}

calibration_fi <- bind_rows(
  make_calibration_data_binary(ukb_bl_imp, "frailty_prevalent", "pred_fi_reference") %>%
    mutate(cohort = "UKB", model = "Reference"),
  make_calibration_data_binary(ukb_bl_imp, "frailty_prevalent", "pred_trim_fi") %>%
    mutate(cohort = "UKB", model = "TRIM-FI"),

  make_calibration_data_binary(esther_bl_imp, "frailty_prevalent", "pred_fi_reference") %>%
    mutate(cohort = "ESTHER", model = "Reference"),
  make_calibration_data_binary(esther_bl_imp, "frailty_prevalent", "pred_trim_fi") %>%
    mutate(cohort = "ESTHER", model = "TRIM-FI")
)

readr::write_csv(
  calibration_fi,
  file.path(dir_tables, "calibration_prevalent_frailty.csv")
)

# ---- 10.11 Decision curve analysis for prevalent frailty ----

run_dca_binary <- function(data, outcome, pred_ref, pred_new, cohort, outcome_label) {
  dat <- data %>%
    transmute(
      outcome = .data[[outcome]],
      Reference = .data[[pred_ref]],
      TRIM = .data[[pred_new]]
    ) %>%
    filter(complete.cases(.))

  dca_fit <- rmda::decision_curve(
    outcome ~ Reference + TRIM,
    data = dat,
    family = binomial(link = "logit"),
    thresholds = seq(0.01, 0.50, by = 0.01),
    confidence.intervals = FALSE,
    study.design = "cohort"
  )

  dca_fit$derived.data %>%
    mutate(
      cohort = cohort,
      outcome_label = outcome_label
    )
}

dca_fi <- bind_rows(
  run_dca_binary(
    data = ukb_bl_imp,
    outcome = "frailty_prevalent",
    pred_ref = "pred_fi_reference",
    pred_new = "pred_trim_fi",
    cohort = "UKB",
    outcome_label = "Prevalent frailty"
  ),
  run_dca_binary(
    data = esther_bl_imp,
    outcome = "frailty_prevalent",
    pred_ref = "pred_fi_reference",
    pred_new = "pred_trim_fi",
    cohort = "ESTHER",
    outcome_label = "Prevalent frailty"
  )
)

readr::write_csv(
  dca_fi,
  file.path(dir_tables, "decision_curve_prevalent_frailty.csv")
)


# ============================================================
# 11. Sensitivity analysis
# ============================================================

# ---- 11.1 Analysis objective ----
# Repeat biomarker selection and prediction analysis using incident
# frailty as the feature-selection outcome in UKB.
#
# Manuscript Supplementary Table 7:
#   Incident-frailty-derived sensitivity panel:
#
#   Proteins:
#     CCL3, FGF23, HGF, IL12B, IL6, LIFR, TGFA, TNFRSF9, VEGFA
#
#   Metabolites:
#     Albumin, Glucose, GlycA, Lactate, S-HDL-PL
#
# Note:
# Supplementary Table 7 text says "6 metabolites", but the footnote and
# Supplementary Table 8 list 5 metabolites. The code follows the listed
# biomarkers.

# ---- 11.2 Sensitivity model variable sets ----

sensitivity_variable_sets <- list(
  Reference = reference_covariates,
  Reference_plus_Proteomics = c(reference_covariates, sensitivity_analysis$proteins),
  Reference_plus_Metabolomics = c(reference_covariates, sensitivity_analysis$metabolites),
  Combined = c(reference_covariates, sensitivity_analysis$omics)
)

# ---- 11.3 Ensure sensitivity biomarkers are available ----

check_required_vars(
  ukb_fp_data,
  sensitivity_variable_sets$Combined,
  "UKB follow-up sensitivity analysis"
)

check_required_vars(
  esther_fp_data,
  sensitivity_variable_sets$Combined,
  "ESTHER follow-up sensitivity analysis"
)

# ---- 11.4 Fit sensitivity TRIM-FP models in UKB ----
# These models correspond to Supplementary Table 7.

sensitivity_fp_models_ukb <- list(
  Reference = fit_ciregic_model(
    data = ukb_fp_data,
    predictors = sensitivity_variable_sets$Reference
  ),

  Reference_plus_Proteomics = fit_ciregic_model(
    data = ukb_fp_data,
    predictors = sensitivity_variable_sets$Reference_plus_Proteomics
  ),

  Reference_plus_Metabolomics = fit_ciregic_model(
    data = ukb_fp_data,
    predictors = sensitivity_variable_sets$Reference_plus_Metabolomics
  ),

  Combined = fit_ciregic_model(
    data = ukb_fp_data,
    predictors = sensitivity_variable_sets$Combined
  )
)

saveRDS(
  sensitivity_fp_models_ukb,
  file.path(dir_models, "sensitivity_incident_selected_fp_models_ukb.rds")
)

# ---- 11.5 Apply sensitivity model coefficients to UKB and ESTHER ----

beta_sens_reference <- stats::coef(sensitivity_fp_models_ukb$Reference)
beta_sens_proteomics <- stats::coef(sensitivity_fp_models_ukb$Reference_plus_Proteomics)
beta_sens_metabolomics <- stats::coef(sensitivity_fp_models_ukb$Reference_plus_Metabolomics)
beta_sens_combined <- stats::coef(sensitivity_fp_models_ukb$Combined)

ukb_fp_data <- ukb_fp_data %>%
  mutate(
    lp_sens_reference = linear_predictor_from_beta(
      data = ukb_fp_data,
      predictors = sensitivity_variable_sets$Reference,
      beta = beta_sens_reference
    ),
    lp_sens_proteomics = linear_predictor_from_beta(
      data = ukb_fp_data,
      predictors = sensitivity_variable_sets$Reference_plus_Proteomics,
      beta = beta_sens_proteomics
    ),
    lp_sens_metabolomics = linear_predictor_from_beta(
      data = ukb_fp_data,
      predictors = sensitivity_variable_sets$Reference_plus_Metabolomics,
      beta = beta_sens_metabolomics
    ),
    lp_sens_combined = linear_predictor_from_beta(
      data = ukb_fp_data,
      predictors = sensitivity_variable_sets$Combined,
      beta = beta_sens_combined
    )
  )

esther_fp_data <- esther_fp_data %>%
  mutate(
    lp_sens_reference = linear_predictor_from_beta(
      data = esther_fp_data,
      predictors = sensitivity_variable_sets$Reference,
      beta = beta_sens_reference
    ),
    lp_sens_proteomics = linear_predictor_from_beta(
      data = esther_fp_data,
      predictors = sensitivity_variable_sets$Reference_plus_Proteomics,
      beta = beta_sens_proteomics
    ),
    lp_sens_metabolomics = linear_predictor_from_beta(
      data = esther_fp_data,
      predictors = sensitivity_variable_sets$Reference_plus_Metabolomics,
      beta = beta_sens_metabolomics
    ),
    lp_sens_combined = linear_predictor_from_beta(
      data = esther_fp_data,
      predictors = sensitivity_variable_sets$Combined,
      beta = beta_sens_combined
    )
  )

# ---- 11.6 Sensitivity C-statistics ----

sensitivity_cstat_results <- bind_rows(
  calc_cstat_incident(ukb_fp_data, "lp_sens_reference", "UKB", "Reference"),
  calc_cstat_incident(ukb_fp_data, "lp_sens_proteomics", "UKB", "Reference + Proteomics"),
  calc_cstat_incident(ukb_fp_data, "lp_sens_metabolomics", "UKB", "Reference + Metabolomics"),
  calc_cstat_incident(ukb_fp_data, "lp_sens_combined", "UKB", "Combined model"),

  calc_cstat_incident(esther_fp_data, "lp_sens_reference", "ESTHER", "Reference"),
  calc_cstat_incident(esther_fp_data, "lp_sens_proteomics", "ESTHER", "Reference + Proteomics"),
  calc_cstat_incident(esther_fp_data, "lp_sens_metabolomics", "ESTHER", "Reference + Metabolomics"),
  calc_cstat_incident(esther_fp_data, "lp_sens_combined", "ESTHER", "Combined model")
)

readr::write_csv(
  sensitivity_cstat_results,
  file.path(dir_tables, "supplementary_table7_sensitivity_cstat.csv")
)

# ---- 11.7 Sensitivity delta C-statistics ----

sensitivity_delta_cstat <- bind_rows(
  bootstrap_delta_cstat(ukb_fp_data, "lp_sens_reference", "lp_sens_proteomics") %>%
    mutate(cohort = "UKB", comparison = "Reference + Proteomics vs Reference"),
  bootstrap_delta_cstat(ukb_fp_data, "lp_sens_reference", "lp_sens_metabolomics") %>%
    mutate(cohort = "UKB", comparison = "Reference + Metabolomics vs Reference"),
  bootstrap_delta_cstat(ukb_fp_data, "lp_sens_reference", "lp_sens_combined") %>%
    mutate(cohort = "UKB", comparison = "Combined model vs Reference"),

  bootstrap_delta_cstat(esther_fp_data, "lp_sens_reference", "lp_sens_proteomics") %>%
    mutate(cohort = "ESTHER", comparison = "Reference + Proteomics vs Reference"),
  bootstrap_delta_cstat(esther_fp_data, "lp_sens_reference", "lp_sens_metabolomics") %>%
    mutate(cohort = "ESTHER", comparison = "Reference + Metabolomics vs Reference"),
  bootstrap_delta_cstat(esther_fp_data, "lp_sens_reference", "lp_sens_combined") %>%
    mutate(cohort = "ESTHER", comparison = "Combined model vs Reference")
)

readr::write_csv(
  sensitivity_delta_cstat,
  file.path(dir_tables, "supplementary_table7_sensitivity_delta_cstat.csv")
)

# ---- 11.8 Save sensitivity datasets with scores ----

readr::write_csv(
  ukb_fp_data,
  file.path(dir_data_processed, "ukb_followup_with_main_and_sensitivity_scores.csv")
)

readr::write_csv(
  esther_fp_data,
  file.path(dir_data_processed, "esther_followup_with_main_and_sensitivity_scores.csv")
)

# ---- 11.9 Sensitivity coefficients ----

sensitivity_coefficients <- purrr::imap_dfr(
  sensitivity_fp_models_ukb,
  extract_ciregic_coefficients
)

readr::write_csv(
  sensitivity_coefficients,
  file.path(dir_tables, "sensitivity_incident_selected_coefficients_ukb.csv")
)
