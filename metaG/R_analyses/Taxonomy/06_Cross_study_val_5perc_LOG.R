## ----setup
list.of.packages <- c("dplyr", "vegan", "BiocManager", "mlr3")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, 
                                          lib = .libPaths()[2],
                                          destdir = .libPaths()[2],
                                          repos="http://cran.us.r-project.org")


list.of.bioc <- c("phyloseq", "SIAMCAT")
new.packages <- list.of.bioc[!(list.of.bioc %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages,
                                          lib = .libPaths()[2],
                                          destdir = .libPaths()[2],)
library(phyloseq)
library(dplyr)
library(SIAMCAT)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (R_analyses dir).n", call.=FALSE)
}


source(paste0(args, "/Scripts/Siamcat_wf.r"))
source(paste0(args, "/Scripts/Cross_study_validation.r"))
set.seed(56987)


## ------------------------------------------------------------------------------
all.tss<-readRDS(paste0(args, "/RDS/all_no0_TSS_noUnassigned.rds"))
all.tss.store<-all.tss 
meta<-readRDS(paste0(args, "/RDS/meta.rds"))


## ------------------------------------------------------------------------------
##############
#  LOG
##############
load(paste0(args, "/RData/metaG_Sinlge_5perc_log.RData"))


## ------------------------------------------------------------------------------
study<-as.character(meta$Study) %>% unique()
spec<-readRDS(paste0(args, "/RDS/spec_5x10.rds")) 
spec<-spec[spec != "unassigned"]
all.tss<-all.tss.store[spec,]
rownames(meta)<-meta$SampleID


print("CSV with lasso")
l_csv_l5<-cross.study.validation(lasso_5x10,
                       microbiome = all.tss,
                       meta=meta,
                       names.studies = unique(meta$Study),
                       name.study.var = "Study",
			is.model.list = T)

print("CSV with enet")
l_csv_e5<-cross.study.validation(enet_5x10,
                       microbiome = all.tss,
                       meta=meta,
                       names.studies = unique(meta$Study),
                       name.study.var = "Study",
			is.model.list = T)

print("CSV with ridge")
l_csv_r5<-cross.study.validation(ridge_5x10,
                       microbiome = all.tss,
                       meta=meta,
                       names.studies = unique(meta$Study),
                       name.study.var = "Study",
			is.model.list = T)

print("CSV with ridge_ll")
l_csv_rl5<-cross.study.validation(ridge_ll_5x10,
                       microbiome = all.tss,
                       meta=meta,
                       names.studies = unique(meta$Study),
                       name.study.var = "Study",
			is.model.list = T)

print("CSV with lasso_ll")
l_csv_ll5<-cross.study.validation(lassoll_5x10,
                       microbiome = all.tss,
                       meta=meta,
                       names.studies = unique(meta$Study),
                       name.study.var = "Study",
			is.model.list = T)

print("CSV with random forest")
l_csv_rf5<-cross.study.validation(rf_5x10,
                       microbiome = all.tss,
                       meta=meta,
                       names.studies = unique(meta$Study),
                       name.study.var = "Study",
			is.model.list = T)


## ------------------------------------------------------------------------------

merge.csv.auc<-function(
  list.csv,
  ML, 
  norm) {
  df<-data.frame()
  for(i in 1:length(list.csv)){
    list.auc<-lapply(list.csv[[i]], function(x) x$summary)
    merged.scv<-do.call(rbind, list.auc)
    df<-rbind.data.frame(df, merged.scv)
  }
  df$ML<-rep(ML, nrow(df))
  df$norm<-rep(norm, nrow(df))

  return(df)
}

csv_log<-rbind(
  merge.csv.auc(l_csv_r5, "ridge", "LOG"),
  merge.csv.auc(l_csv_e5, "enet", "LOG"),
  merge.csv.auc(l_csv_l5, "lasso", "LOG"),
  merge.csv.auc(l_csv_rl5, "ridge_ll", "LOG"),
  merge.csv.auc(l_csv_ll5, "lasso_ll", "LOG"),
  merge.csv.auc(l_csv_rf5, "random_forest", "LOG")

  )

saveRDS(csv_log, paste0(args, "/RDS/metaG_CSV_5perc_Log_df_auc.rds"))


## ------------------------------------------------------------------------------
print("Save wsp")
save.image(paste0(args, "/RData/metaG_Cross_validation_5perc_LOG.Rdata"))

