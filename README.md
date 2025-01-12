### Overview

These scripts have been used to conduct the data-analysis in the publication:


Machine learning-based meta-analysis reveals gut microbiome alterations associated with Parkinson’s disease
Stefano Romano, Jakob Wirbel, Rebecca Ansorge, Christian Schudoma, Quinten Raymond Ducarmon, Arjan Narbad, Georg Zeller
bioRxiv 2023.12.05.569565; doi: https://doi.org/10.1101/2023.12.05.569565 

This code was tailored to the analyses conducted in the above manuscript, and it is not intended to be used in other contexts.
The work-flow consists of a series of numbered R scripts that need to be executed one after the other. There are additional R scripts (e.g. AA_, AB_, AC_) that were used to create the final figures and tables.

If you find this code useful or use any of the approaches below in your work, please cite:

Machine learning-based meta-analysis reveals gut microbiome alterations associated with Parkinson’s disease
Stefano Romano, Jakob Wirbel, Rebecca Ansorge, Christian Schudoma, Quinten Raymond Ducarmon, Arjan Narbad, Georg Zeller
bioRxiv 2023.12.05.569565; doi: https://doi.org/10.1101/2023.12.05.569565 

### Taxonomic and functional profiling

#### Taxonomic profiling of 16S amplicon sequences
16S amplicon seqeunces were profiled using DADA2 following the authors guidelines accessible at https://benjjneb.github.io/dada2/tutorial.html
Each dataset was profiled indipendently. In addition, samples sequenced on different runs were profiled independently to allow a run-specific estimation of the sequencing error rates.

#### Taxonomic profiling of shotgun metagenomic sequences
Taxonomic profiling of the shotgun metagenomic sequences were performed using mOTUs_v3 https://github.com/motu-tool/mOTUs.
Each dataset was profile indipendently using the command:

e.g. `motus profile -s [NAME_INPUT].fastq.gz -o [NAME_OUTPUT].v3.motus -n [NAME_INPUT] -g 2 -l [LENGHT e.g. 50] -c`

Samples profiles were then merged using:

e.g. `motus merge -i [NAME_INPUT] > [NAME_OUTPUT]_all_samples.v3.g2_l75.motus`

#### Functional profiling of shotgun metagenomic sequences
Quality and host-filtered reads were aligned to the reduced GMGC (downlodabel from https://github.com/zellerlab/cayman).
Resulting bam files were then processed using gffquant (https://github.com/cschu/gff_quantifier) to obtain final abundance tables.


### Repository strucutre

The structure of the repository is as follow:
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
│   │   ├── 14_Coefficients.Rmd
│   │   ├── 15_DA_16S.Rmd
│   │   ├── 16_DA_Figure.Rmd
│   │   ├── 17_Sex_Age_cov.rmd
│   │   ├── 18_extract_AUCs_Iters.Rmd
│   │   ├── 19_Jo.match16S.SMG.rmd
│   │   ├── AA_Figures.Rmd
│   │   ├── AB_CSV_Figures.Rmd
│   │   └── AC_metaG_16S.Rmd
│   └── R_analyses_AD_MS
│       ├── 01_Filtration.Rmd
│       ├── 02_CDV.Rmd
│       ├── 03_CDV_LOSO.Rmd
│       └── 04_Figure_CDV.Rmd
├── LICENSE
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
│           ├── 10_Combinatorics_LOSO.Rmd
│           ├── 11_Combinatorics_LOSO_figs.Rmd
│           ├── 12_extract_AUCs_Iters.Rmd
│           ├── 13_Coeff_vs_DA.Rmd
│           ├── 14_LOSO_feat.sel_blocking.Rmd
│           ├── 17_LOSO_feat.sel_figs.Rmd
│           ├── AA_Figures.Rmd
│           └── AB_CSV_metaG.Rmd
└── Scripts
    ├── compare.groups.r
    ├── Cross_study_validation.r
    ├── DA_function_gen_odd.R
    ├── LOSO.r
    ├── Model_coefficient_extraction.r
    ├── Prevalence_functions.r
    └── Siamcat_wf.r
```

### Packages installation

Installation of the R packages used in the analyses is reported in the respective R scripts. In general this has been done using either `install.packages("package_name")` or `BiocManager::install("package_name")`

### Usage instruction

The workflow for each profile is as follows:

   * filter the datasets
   * perform beta diversity analyses (this has been done only for the taxonomic profiles derived from the shotgun metagenomics and 16S amplicon data)
   * rarefy and batch correct the data (only for 16S amplicon data)
   * run machine learning models for each dataset
   * perform cross study validation (CSV)
   * perform leave-one-study out validation (LOSO)
   * perform combinatorics experiment using LOSO models (only for metagenomics data)
   * perform feature selection and new LOSO validation (only for metagenomics data)
   * perform cross disease prediction (CDV; only for 16S amplicon data)
   * extract model coefficients (for the taxonomy profiles the coefficients have been used to build ordinations)
   * perform differential abundance (DA) analyses
   * perform confounders analyses
   * perform enrichment analysis for KO/KEGG pathways 
   * match sample IDs between 16S and metagenomics data (possible only in one dataset) and run new ML models matching training and test splits between data types
   * combine results in figures and tables

ML models built in the KO data were run in two different scripts (03a, 03b; 05a, 05b) to accommodate memory restrictions during computation.

### Package versions

All data analyses were performed in R v_4.2 and using Bioconductor v_3.15. We used the following R packages: `omixerRpm v_0.3.3`, `phyloseq v_1.40`, `vegan v_2.6.4`, `genodds v_1.1.2`, `meta v_6.2.1`, `leaps v_3.1`, `SIAMCAT v_2.0`, `SIAMCAT v_2.10`, `stats v_4.2.3`, `rtk v_0.2.6.1`, `MMUPHin v_1.10.3`, `bapred v_1.1`, `nlme v_3.1.162`, `emmeans v_1.8.5`, `performance v_0.11`, `ggplot2 v_3.4.4`, `mlr3 v_0.15`, `pROC v_1.18.0`, `mlr3extralearners v_0.6`, `coin v_1.4.2`, `ggnewscale v_0.4.9`, `cowplot v_1.1.1`, `BiocManager v_1.30.20`, `clusterProfiler v_4.4.4` `rstatix v_0.7.2`, `ggh4x v_0.2.5`, `RColorBrewer v_1.1.3`, `microbiome v_1.18.0`



