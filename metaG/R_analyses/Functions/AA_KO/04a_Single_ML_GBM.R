## ----setup, include=FALSE-----------------------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------------------------------------------
all.tss<-readRDS(paste0(args, "/RDS/GBM.rds"))
meta<-readRDS(paste0(args, "/RDS/meta.rds"))


## -----------------------------------------------------------------------------------------------------------------
study<-as.character(meta$Study) %>% unique()
rownames(meta)<-meta$SampleID
all.tss[is.na(all.tss)]<-0
min<-min(unique(as.vector(as.matrix(all.tss)))[unique(as.vector(as.matrix(all.tss))) != 0])/100

print("Running Lasso")
lasso_5x10<-siamcat.wf.lists(vector.names = study,
                       numCores = 10, df = all.tss, meta = meta, name.feat = "PD",  filt.method = "pass", 
                       cutoff = 0.30, log.n0 = min, nsplit = 1, n.fold = 10, 
                       trafo = "log.std",  ml.method = "lasso", 
                       raw.auc = F, plots = F,
                       max.show = 100, name.study.var = "Study", reps = 1, 
                       perform.selection = F,
                       param.selection = list(no_features = 600, method = "AUC", direction="absolute"))
auc<-lasso_5x10@auc
lassoo.df<-data.frame(study = names(auc), prev = rep("5% in 2 std", 
                    length(lasso_5x10)), ML = rep("lasso", length(auc)),
                    AUC = unlist(lapply(auc, function(x) x$auc %>% as.double())),
                    norm = rep("log", length(auc)))

print("Running lasso_ll")
lassoll_5x10<-siamcat.wf.lists(vector.names = study,
                       numCores = 10, df = all.tss, meta = meta, name.feat = "PD",  filt.method = "pass", 
                       cutoff = 0.30, log.n0 = min, nsplit = 1, n.fold = 10, 
                       trafo = "log.std",  ml.method = "lasso_ll", 
                       raw.auc = F, plots = F,
                       max.show = 100, name.study.var = "Study", reps = 1,
                       perform.selection = F, 
                       param.selection = list(no_features = 600, method = "AUC", direction='absolute'))
auc<-lassoll_5x10@auc
lassoll.df<-data.frame(study = names(auc), prev = rep("5% in 2 std", 
                    length(lasso_5x10)), ML = rep("lasso_ll", length(auc)),
                    AUC = unlist(lapply(auc, function(x) x$auc %>% as.double())),
                    norm = rep("log", length(auc)))

print("Running ridge")
ridge_5x10<-siamcat.wf.lists(vector.names = study,
                       numCores = 10, df = all.tss, meta = meta, name.feat = "PD",  filt.method = "pass", 
                       cutoff = 0.30, log.n0 = min, nsplit = 1, n.fold = 10, 
                       trafo = "log.std",  ml.method = "ridge", 
                       raw.auc = F, plots = F,
                       max.show = 100, name.study.var = "Study", reps = 1,
                       perform.selection = F, 
                       param.selection = list(no_features = 600, method = "AUC", direction='absolute'))
auc<-ridge_5x10@auc
ridge.df<-data.frame(study = names(auc), prev = rep("5% in 2 std", 
                    length(auc)), ML = rep("ridge", length(auc)),
                    AUC = unlist(lapply(auc, function(x) x$auc %>% as.double())),
                    norm = rep("log", length(auc)))


print("Running ridge_ll")
ridge_ll_5x10<-siamcat.wf.lists(vector.names = study,
                       numCores = 10, df = all.tss, meta = meta, name.feat = "PD",  filt.method = "pass", 
                       cutoff = 0.30, log.n0 = min, nsplit = 1, n.fold = 10, 
                       trafo = "log.std",  ml.method = "ridge_ll", 
                       raw.auc = F, plots = F,
                       max.show = 100, name.study.var = "Study", reps = 1,
                       perform.selection = F, 
                       param.selection = list(no_features = 600, method = "AUC", direction='absolute'))
auc<-ridge_ll_5x10@auc
ridgell.df<-data.frame(study = names(auc), prev = rep("5% in 2 std", 
                    length(auc)), ML = rep("ridge_ll", length(auc)),
                    AUC = unlist(lapply(auc, function(x) x$auc %>% as.double())),
                    norm = rep("log", length(auc)))


print("Running enet")
enet_5x10<-siamcat.wf.lists(vector.names = study,
                       numCores = 10, df = all.tss, meta = meta, name.feat = "PD",  filt.method = "pass", 
                       cutoff = 0.30, log.n0 = min, nsplit = 1, n.fold = 10, 
                       trafo = "log.std",  ml.method = "enet", 
                       raw.auc = F, plots = F,
                       max.show = 100, name.study.var = "Study",
                       param.set=list('alpha'=0.5), reps = 1,
                       perform.selection = F, 
                       param.selection = list(no_features = 600, method = "AUC", direction='absolute'))
auc<-enet_5x10@auc
enet.df<-data.frame(study = names(auc), prev = rep("5% in 2 std", 
                    length(auc)), ML = rep("enet", length(auc)),
                    AUC = unlist(lapply(auc, function(x) x$auc %>% as.double())),
                    norm = rep("log", length(auc)))


print("Running randomForest")
rf_5x10<-siamcat.wf.lists(vector.names = study,
                       numCores = 10, df = all.tss, meta = meta, name.feat = "PD",  filt.method = "pass", 
                       cutoff = 0.30, log.n0 = min, nsplit = 1, n.fold = 10, 
                       trafo = "log.std",  ml.method = "randomForest", 
                       raw.auc = F, plots = F,
                       max.show = 100, name.study.var = "Study", reps = 1,
                       perform.selection = F, 
                       param.selection = list(no_features = 600, method = "AUC", direction='absolute'))

auc<-rf_5x10@auc
rf.df<-data.frame(study = names(auc), prev = rep("5% in 2 std", 
                    length(auc)), ML = rep("random_forest", length(auc)),
                   AUC = unlist(lapply(auc, function(x) x$auc %>% as.double())),
                    norm = rep("log", length(auc)))

print("creating asn saving AUC-df")
df<-rbind.data.frame(lassoo.df, enet.df, ridge.df, ridgell.df, lassoll.df, rf.df)
saveRDS(df, paste0(args, "/RDS/GBM_single_5per_log.rds"))


## -----------------------------------------------------------------------------------------------------------------
print("Saving wsp")
save.image(paste0(args, "/RData/GBM_Sinlge_5perc_log.RData"))

