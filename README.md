## PD_ML_meta


The structure of the repos is as follow:
```
├── 16S
│   ├── R_analyses
│   │   ├── 01_Filtration.Rmd
│   │   ├── 02_Rarefy.Rmd
│   │   ├── 03_Beta_16S.Rmd
│   │   ├── 04a_ML_single_study_LOG_5perc.R
│   │   ├── 04b_ML_single_study_LOG_10perc.R
│   │   ├── 04c_ML_single_study_LOG_20perc.R
│   │   ├── 04d_ML_single_study_LOG_30perc.R
│   │   ├── 06a_ML_single_study_CLR_5perc.R
│   │   ├── 06b_ML_single_study_CLR_10perc.R
│   │   ├── 06c_ML_single_study_CLR_20perc.R
│   │   ├── 06d_ML_single_study_CLR_30perc.R
│   │   ├── 07a_ML_single_study_LOG_5perc_rar2000.R
│   │   ├── 08a_Cross_study_val_5perc_LOG.R
│   │   ├── 08b_Cross_study_val_Rarefied_5perc_LOG.R
│   │   ├── 10_Correct_batch.Rmd
│   │   ├── 11a_Batch_MMuphin.R
│   │   ├── 11b_Batch_Mean.R
│   │   ├── 11c_Batch_RatioA.R
│   │   ├── 11d_Batch_RatioG.R
│   │   ├── 12a_MMUphin_CSV.R
│   │   ├── 12b_Mean_CSV.R
│   │   ├── 12c_RatioA_CSV.R
│   │   ├── 12d_RatioG_CSV.R
│   │   ├── 13_ML_loso_ridge.R
│   │   ├── 14_Coefficiens.Rmd
│   │   ├── 15_DA_16S.Rmd
│   │   ├── 16_DA_Figure.Rmd
│   │   ├── 17_Sex_Age_cov.rmd
│   │   ├── AA_Figures.Rmd
│   │   ├── AB_CSV_Figures.Rmd
│   │   └── AC_16S_metaG_Fig2.Rmd
│   └── R_analyses_AD_MS
│       ├── 01_Filtration.Rmd
│       ├── 02_CDV.Rmd
│       ├── 03_CDV_LOSO.Rmd
│       └── 04_Figure_CDV.Rmd
├── metaG
│   └── R_analyses
│       ├── Functions
│       │   ├── AA_KO
│       │   │   ├── 01_Filter_datasets.Rmd
│       │   │   ├── 03a_ML_single_study_LOG_5perc.R
│       │   │   ├── 03b_Wallen_ML_single_study_LOG_5perc.R
│       │   │   ├── 04a_Single_ML_GBM.R
│       │   │   ├── 04b_Single_ML_GMM.R
│       │   │   ├── 05a_Cross_study_val_5perc_LOG.R
│       │   │   ├── 05b_KO_CSV_Wallen.Rmd
│       │   │   ├── 05c_CSV_GMM_5perc_LOG.R
│       │   │   ├── 05d_CSV_GBM_5perc_LOG.R
│       │   │   ├── 06a_ML_loso_GMM.R
│       │   │   ├── 06b_ML_loso_GBM.R
│       │   │   ├── 06_ML_loso.R
│       │   │   ├── 07_Coefficients.Rmd
│       │   │   ├── 08_DA_KO.Rmd
│       │   │   ├── 09a_DA_figure_GBM.Rmd
│       │   │   ├── 09b_DA_figure_GMM.Rmd
│       │   │   ├── 10_DA_table_KO.Rmd
│       │   │   ├── 11_Enrichment_KO.Rmd
│       │   │   ├── AA_Figures.Rmd
│       │   │   ├── AB_Violin.Rmd
│       │   │   └── AC_Network.Rmd
│       │   ├── AB_KOmod
│       │   │   ├── 01_Filtration.Rmd
│       │   │   ├── 02_ML_single_study_LOG_5perc.R
│       │   │   ├── 03_Cross_study_val_5perc_LOG.R
│       │   │   ├── 04_ML_loso.R
│       │   │   ├── 05_Coefficients.Rmd
│       │   │   ├── 06_DA_KOmod.Rmd
│       │   │   ├── 07_DA_tab_Kmod.Rmd
│       │   │   └── AA_Figures.Rmd
│       │   └── AC_KOpath
│       │       ├── 01_Filtration.Rmd
│       │       ├── 03_ML_single_study_LOG_5perc.R
│       │       ├── 04_Cross_study_val_5perc_LOG.R
│       │       ├── 05_ML_loso.R
│       │       ├── 06_Coefficients.Rmd
│       │       ├── 07_DA_KOpath.Rmd
│       │       ├── 08_DA_figure.Rmd
│       │       ├── 09_Cov_drugs.rmd
│       │       ├── 10_Cov_drugs_table.rmd
│       │       ├── 11_Cov_PD_drugs.rmd
│       │       ├── 12_Cov_PD_drugs_table.rmd
│       │       ├── 13_Cov_sex_age.rmd
│       │       ├── 14_Cov_figs.rmd
│       │       └── AA_Figures.Rmd
│       └── Taxonomy
│           ├── 01_Filtration.Rmd
│           ├── 02_Beta.Rmd
│           ├── 03_DA_metaG.Rmd
│           ├── 05a_ML_single_study_CLR_5perc.R
│           ├── 05a_ML_single_study_LOG_5perc.R
│           ├── 05b_ML_single_study_CLR_10perc.R
│           ├── 05b_ML_single_study_LOG_10perc.R
│           ├── 05c_ML_single_study_CLR_20perc.R
│           ├── 05c_ML_single_study_LOG_20perc.R
│           ├── 05d_ML_single_study_CLR_30perc.R
│           ├── 05d_ML_single_study_LOG_30perc.R
│           ├── 06_Cross_study_val_5perc_LOG.R
│           ├── 07_ML_loso_v1.R
│           ├── 08_Coefficients.Rmd
│           ├── 09_DA_figure.Rmd
│           ├── AA_Figures.Rmd
│           ├── AB_CSV_metaG.Rmd
│           └── AC_metaG_16S.Rmd
└── Scripts
    ├── Cross_study_validation.r
    ├── DA_function_gen_odd.R
    ├── Model_coefficient_extraction.r
    ├── Pevalence_functions.r
    └── Siamcat_wf.r
```
The workflow for each profile is as follows:

   * filter the datasets
   * perform beta diversity analyses (this has been done only for the taxonomic profiles derived from the shotgun metagenomics and 16S amplicon data)
   * rarefy and batch correct the data (only for 16S amplicon data)
   * run machine learning models for each dataset
   * perform cross study validation (CSV)
   * perform leave-one-study out validation (LOSO)
   * perform cross disease prediction (CDV; only for 16S amplicon data)
   * extract model coefficients (for the taxonomy profiles the coefficients have been used to build ordinations)
   * perform differential abundance (DA) analyses
   * perform confounders analyses
   * combine results in figures and tables

Analyses for shot-gun metagenomes (metaG) were run first. This produced some intermediate files that were then used while processing 16S amplicon data for producing the final plots. ML models built in the KO data were run in two different scripts (03a, 03b; 05a, 05b) to accommodate memory restrictions during computation.

