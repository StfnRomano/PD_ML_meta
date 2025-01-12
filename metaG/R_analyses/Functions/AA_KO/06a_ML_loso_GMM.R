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
require(mlr3extralearners)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (R_analyses dir).n", call.=FALSE)
}


source(paste0(args, "/scripts/Siamcat_wf.r"))
set.seed(56987)


## -----------------------------------------------------------------------------
all.tss<-readRDS(paste0(args, "/RDS/GMM.rds"))
all.tss.store<-all.tss 
meta<-readRDS(paste0(args, "/RDS/meta.rds"))


## -----------------------------------------------------------------------------
study<-as.character(meta$Study) %>% unique()
min<-min(unique(as.vector(as.matrix(all.tss)))[unique(as.vector(as.matrix(all.tss))) != 0])/100
rownames(meta)<-meta$SampleID
meta<-meta[order(meta$SampleID),]


loso_ridge<-loso.validation( vector.names = study, 
                           df = all.tss, meta = meta, name.feat = "PD",
                           filt.method = "pass", cutoff = 1e-07, trafo = "log.std", 
                           log.n0 = min, n.fold = 10, 
                           10, ml.method = "ridge", 
                           name.study = NULL, dir = NULL, 
                           max.show = 50, feat.type = "filtered", sd.add = 0.1, raw.auc = F, plots = F,
                           siamcat.holdout = NULL, param.set=NULL,
                           study.label = "Study",
                           case = "PD",
                           n.p=2,
                           norm.margin = 1,
			   perform.selection = FALSE, 
                       	   param.selection = list(no_features = 600, method = "AUC", direction="absolute")) 
saveRDS(loso_ridge, paste0(args, "/RDS/GMM_loso_ridge.rds"))


