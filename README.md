# TRIM Frailty Analysis Code

This repository contains the main R analysis code for the TRIM frailty manuscript.

The code implements the derivation and external validation of two omics-based frailty models:

- **TRIM-FI**: Traditional Risk factor, Inflammation, and Metabolism based Frailty Identification model
- **TRIM-FP**: Traditional Risk factor, Inflammation, and Metabolism based Frailty Prediction model

The main analysis uses **UK Biobank** as the derivation cohort and **ESTHER** as the external validation cohort.

---

## Repository overview

The main analysis script covers:

- Project setup
- Biomarker panels and conventional covariates
- Data loading
- Frailty index construction
- Omics preprocessing
- Missing data imputation
- Prevalent frailty LASSO analysis
- TRIM-FI model fitting and validation
- TRIM-FP model fitting and validation
- External validation in ESTHER
- Model performance evaluation
- Sensitivity analysis

The script focuses on data preparation, model fitting, external validation, performance evaluation, and sensitivity analysis. Manuscript tables and figures can be generated separately from the saved outputs.

---

## Contact

For questions about the analysis code, please contact the first author:

```text
penglei960825@gmail.com
```
