## ----setup, include=FALSE------------------------------------------------------------------------
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

# load my function to run siamcat
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (R_analyses dir).n", call.=FALSE)
}


source(paste0(args, "/Scripts/Siamcat_wf.r"))
set.seed(56987)


## ------------------------------------------------------------------------------------------------
all.tss<-readRDS(paste0(args, "/RDS/all_no0_TSS_KO.rds"))
all.tss.store<-all.tss 
meta<-readRDS(paste0(args, "/RDS/meta.rds"))

# just select Wallen here.
st<-unique(meta$Study)[7]
meta<-subset(meta, Study %in% st)
all.tss.store<-all.tss[, which(names(all.tss) %in% meta$SampleID)]


## ------------------------------------------------------------------------------------------------
study<-as.character(meta$Study) %>% unique()
spec<-readRDS(paste0(args, "/RDS/spec_5x10.rds"))
all.tss<-all.tss.store[spec,]
rownames(meta)<-meta$SampleID

min<-min(unique(as.vector(as.matrix(all.tss)))[unique(as.vector(as.matrix(all.tss))) != 0])/100

# I do not colelct the AUC because this are single models
print("Running Lasso")
lasso_5x10<-siamcat.wf.lists(vector.names = study,
                       numCores = 0, df = all.tss, meta = meta, name.feat = "PD",  filt.method = "pass", 
                       cutoff = 0.30, log.n0 = min, nsplit = 10, n.fold = 10, 
                       trafo = "log.std",  ml.method = "lasso", 
                       raw.auc = F, plots = F,
                       max.show = 100, name.study.var = "Study", reps = 0, 
                       perform.selection = TRUE,
                       param.selection = list(no_features = 600, method = "AUC", direction="absolute"))

print("Running lasso_ll")
lassoll_5x10<-siamcat.wf.lists(vector.names = study,
                       numCores = 0, df = all.tss, meta = meta, name.feat = "PD",  filt.method = "pass", 
                       cutoff = 0.30, log.n0 = min, nsplit = 10, n.fold = 10, 
                       trafo = "log.std",  ml.method = "lasso_ll", 
                       raw.auc = F, plots = F,
                       max.show = 100, name.study.var = "Study", reps = 0,
                       perform.selection = TRUE, 
                       param.selection = list(no_features = 600, method = "AUC", direction='absolute'))

print("Running ridge")
ridge_5x10<-siamcat.wf.lists(vector.names = study,
                       numCores = 0, df = all.tss, meta = meta, name.feat = "PD",  filt.method = "pass", 
                       cutoff = 0.30, log.n0 = min, nsplit = 10, n.fold = 10, 
                       trafo = "log.std",  ml.method = "ridge", 
                       raw.auc = F, plots = F,
                       max.show = 100, name.study.var = "Study", reps = 0,
                       perform.selection = TRUE, 
                       param.selection = list(no_features = 600, method = "AUC", direction='absolute'))


print("Running ridge_ll")
ridge_ll_5x10<-siamcat.wf.lists(vector.names = study,
                       numCores = 0, df = all.tss, meta = meta, name.feat = "PD",  filt.method = "pass", 
                       cutoff = 0.30, log.n0 = min, nsplit = 10, n.fold = 10, 
                       trafo = "log.std",  ml.method = "ridge_ll", 
                       raw.auc = F, plots = F,
                       max.show = 100, name.study.var = "Study", reps = 0,
                       perform.selection = TRUE, 
                       param.selection = list(no_features = 600, method = "AUC", direction='absolute'))


print("Running enet")
enet_5x10<-siamcat.wf.lists(vector.names = study,
                       numCores = 0, df = all.tss, meta = meta, name.feat = "PD",  filt.method = "pass", 
                       cutoff = 0.30, log.n0 = min, nsplit = 10, n.fold = 10, 
                       trafo = "log.std",  ml.method = "enet", 
                       raw.auc = F, plots = F,
                       max.show = 100, name.study.var = "Study",
                       param.set=list('alpha'=0.5), reps = 0,
                       perform.selection = TRUE, 
                       param.selection = list(no_features = 600, method = "AUC", direction='absolute'))


print("Running randomForest")
rf_5x10<-siamcat.wf.lists(vector.names = study,
                       numCores = 0, df = all.tss, meta = meta, name.feat = "PD",  filt.method = "pass", 
                       cutoff = 0.30, log.n0 = min, nsplit = 10, n.fold = 10, 
                       trafo = "log.std",  ml.method = "randomForest", 
                       raw.auc = F, plots = F,
                       max.show = 100, name.study.var = "Study", reps = 0,
                       perform.selection = TRUE, 
                       param.selection = list(no_features = 600, method = "AUC", direction='absolute'))


## ------------------------------------------------------------------------------------------------
print("Saving wsp")
save.image(paste0(args, "/RData/metaG_Sinlge_5perc_log_Wallen.RData"))

