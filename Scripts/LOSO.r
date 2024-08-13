
loso.combinatorics<-function( list.combination, # list of study combinations
                              test.name, # names of study to use as test
                              df, meta,
                              sd.min.q = 0.1,
                              log.n0,
                              study.label = "Study",
                              make.names = F, smg =T) {
                                #change names of motus
                                if(smg == T){
                                  mot<-rownames(df)
                                  rownames(df)<-gsub("\\]", "", gsub(".* \\[", "", mot))
                                }
                                
                                
                                # function for finding best lambda - from Siamcat
                                get.best.glmnet.lambda <- function(model, measure, min.nonzero = 5, task){
                                  idx <- which(model$model$nzero >= min.nonzero)
                                  new.model <- model$clone(deep = TRUE)
                                  perf <- vapply(idx, FUN=function(x){
                                    new.model$param_set$values$s <- new.model$model$lambda[x]
                                    pred <- new.model$predict(task)
                                    pred$score(msr(measure))
                                  }, FUN.VALUE = double(1))
                                  if (msr(measure)$minimize){
                                    f.idx <- which(perf==min(perf))[1]
                                  } else {
                                    f.idx <- which(perf==max(perf))[1]
                                  }
                                  new.model$param_set$values$s <- new.model$model$lambda[f.idx]
                                  return(new.model)
                                }
                                # normalize the data [as in SIAMCAT]
                                message("---> normilize features with log.std to mimic Siamcat behavior")
                                feat.log <- log10(df + log.n0)
                                m <- rowMeans(feat.log)
                                s <- apply(feat.log, MARGIN = 1, FUN = sd)
                                q <- quantile(s, sd.min.q, names = FALSE)
                                stopifnot(q > 0)
                                feat.norm <- (feat.log - m) / (s + q)
                                
                                
                                l_smct<-vector(mode = "list", length = length(list.combination)) # store for each study a list of models
                                l_smct.test<-vector(mode = "list", length = length(list.combination)) # store for the test siamcat objects
                                l_auroc<-vector(mode = "list", length = length(list.combination))
                                names(l_smct)<-names(list.combination)
                                names(l_smct.test)<-names(list.combination)
                                names(l_auroc)<-names(list.combination)
                                
                                l.df_loso_cv<-vector(mode = "list", length = length(list.combination))
                                #######################################
                                # Define the test object
                                # (it is always the same)
                                ###############################
                                message("Select data fors testing set ", test.name)
                                # make testing set
                                # change names of motus
                                if(make.names == T){
                                  row.names(df) <- make.names(rownames(df))
                                } 
                                
                                tmp.test<-subset(meta,  meta[, study.label]  == test.name)
                                tmp.test.feat<-as.data.frame(t(df[, rownames(tmp.test)]))
                                
                                stopifnot(rownames(tmp.test) == tmp.test$SampleID)
                                tmp.test.feat$PD<-tmp.test$PD %>%
                                  as.factor()
                                
                                # loop over subsets
                                for(i in 1:length(list.combination)){
                                  # extract the subsets for iter i
                                  comb.trains<-as.data.frame(list.combination[i])
                                  # define df and l.smcts for iter i 
                                  df_loso_cv<-data.frame(study.train=vector(length = ncol(comb.trains)), 
                                                         AUC = vector(length = ncol(comb.trains)), 
                                                         n.samples = vector(length = ncol(comb.trains)), 
                                                         n.pd =vector(length = ncol(comb.trains)), 
                                                         n.hc =vector(length = ncol(comb.trains)))
                                  
                                  l_smct.tmp<-vector(mode = "list", length = ncol(comb.trains)) # store for each study a list of models
                                  names(l_smct.tmp)<-names(comb.trains)
                                  l_smct.test.tmp<-vector(mode = "list", length = ncol(comb.trains)) # store for the test siamcat objects
                                  names(l_smct.test.tmp)<-names(comb.trains)
                                  
                                  lc<-names(list.combination)[i]
                                  message("Training of the sets ", lc)
                                  l.tmp<-vector(mode = "list", length = ncol(comb.trains))
                                  for(cmt in 1:ncol(comb.trains)){
                                    cm.n<-names(comb.trains)
                                    message("---> Training on ", cm.n[cmt], ", including studies: ", comb.trains[,cmt])
                                    vector.names<-comb.trains[,cmt]
                                    # make training set
                                    tmp<-subset(meta, meta[, study.label] %in% vector.names)
                                    tmp.feat<-df[, rownames(tmp)] 
                                    # fill in info on sample size in df
                                    df_loso_cv[cmt,1]<-paste(vector.names, collapse = ", ")
                                    df_loso_cv[cmt,3]<-nrow(tmp)
                                    df_loso_cv[cmt,4]<-nrow(subset(tmp, PD == "PD"))
                                    df_loso_cv[cmt,5]<-nrow(subset(tmp, PD == "HC"))
                                    # create training sciamcat object
                                    message("---> define the mlr3 objects")
                                    # define lerner to be class lasso/enet/ridge
                                    learner <- lrn("classif.cv_glmnet")
                                    learner$predict_type <- 'prob'
                                    learner$param_set$values$alpha <- 0 # for ridge
                                    tmp.feat<-as.data.frame(t(tmp.feat))
                                    stopifnot(rownames(tmp.feat) == tmp$SampleID)
                                    tmp.feat$PD<-tmp$PD %>% as.factor()
                                    # define the task
                                    task.train <- TaskClassif$new(id='classif', backend=tmp.feat, target = "PD", positive = "PD")
                                    # use siamcat approach to find best lambda
                                    model <- learner$train(task = task.train)
                                    message("---> select lambda to mimic Siamcat behavior")
                                    model<-get.best.glmnet.lambda(model, measure = "classif.auc", task = task.train)
                                    l_smct[[cmt]]<-model
                                    # make external prediction
                                    # define test task
                                    message("---> defining test task")
                                    test.task <- TaskClassif$new(id='classif', backend=tmp.test.feat,
                                                                 target='PD', positive = "PD")
                                    
                                    message("---> make prediction on holdout data: ", test.name)
                                    pdata <- model$predict(task=test.task)
                                    roc<-pROC::roc(response=tmp.test$PD, predictor = pdata$prob[,1], direction = "<")
                                    # Running  model with holdout
                                    l_smct.test.tmp[[cmt]]<-pdata
                                    
                                    # extract auc
                                    df_loso_cv[cmt,2]<-roc$auc %>% as.double()
                                    
                                  }
                                  l_smct[[i]]<-l_smct.tmp
                                  l_smct.test[[i]]<-l_smct.test.tmp
                                  l.df_loso_cv[[i]]<-df_loso_cv
                                }
                                
                                setClass("siamcat.wf.out", slots=list(siamcat.train="list", siamca.test = "list", auc="list"))
                                out<-new("siamcat.wf.out", siamcat.train = l_smct, siamca.test = l_smct.test, auc = l.df_loso_cv)
                                return(out)
                              }




feats.20w.sel.block.ns<-function(meta,
                                 microbiome,
                                 study.vector = study) {
  require(coin)
  meta$PD<-as.character(meta$PD)
  # create lists
  l.feats<-vector(mode = "list", length = 7)
  l.df<-vector(mode = "list", length = 7)
  
  # name lists
  names(l.feats)<-study.vector
  names(l.df)<-study.vector
  
  message("Check whether the samples are properly formatted")
  stopifnot(all(colnames(microbiome) == meta$SampleID)) # remember taxa are rows
  
  for(i in 1:length(study.vector)){
    
    message("Running selection on LOSO tested on ", study[i])
    message("---> finding the sample ids for train and test set")
    tr<-study[-i]
    te<-study[i]
    train.id<-which(meta$Study %in% tr)
    test.id<-which(meta$Study %in% te)
    message("---> define train set for doing feats selection using Wilcoxon")
    
    meta.train<-meta[train.id,]
    df.train<-microbiome[,which(names(microbiome) %in% meta.train$SampleID)]
    df.train<-as.data.frame(t(df.train))
    stopifnot(rownames(df.train) == meta.train$SampleID)
    message("Size of the meta.train is ", (dim(meta.train)))
    message("Size of the meta.train is ", (dim(df.train)))
    message("The left out study is ", study[i], ", and the training is based on: ")
    print(unique(meta.train$Study))
    
    # transpose df for the wilcoxon
    df.da<-data.frame(w = vector(length = ncol(df.train)),
                      pvalue = vector(length = ncol(df.train)),
                      motus = vector(length = ncol(df.train)))
    message("---> perform Wilcoxon test to select top 20 feats blocking by study")
    for(f in 1:ncol(df.train)){
      tmp<-data.frame(feats = df.train[,f],
                      pd = meta.train$PD %>% as.factor(),
                      study = meta.train$Study %>% as.factor())
      
      w<-coin::wilcox_test(feats ~ pd | study, data = tmp,)
      df.da$w[f]<-coin::statistic(w)
      df.da$pvalue[f]<-coin::pvalue(w)
      df.da$motus[f]<-names(df.train)[f]    
    }
    # sort effect size
    df.da$abs<-abs(df.da$w)
    df.da$padjust<-p.adjust(df.da$pvalue, method = "fdr")
    df.da.sign<-subset(df.da, padjust < 0.05)
    feats.tmp<-head(df.da.sign[order(df.da.sign$abs, decreasing = T),], n = 20)$motus
    
    l.df[[i]]<-df.da
    l.feats[[i]]<-feats.tmp
    
  }
  message("Create S4 object and closing")
  
  setClass("loso.w20",representation(l.df="list",l.f="list"))
  out<-new("loso.w20",l.df=l.df,l.f=l.feats)
  return(out)
}
