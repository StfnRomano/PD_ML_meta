#################################################################################
#
# This is a function to run all SIAMCAT steps 
# At the end, the scripts gives a list including SIAMCAT object + auroc
# 
#
#################################################################################



siamcat.wf<-function(df, meta, name.feat, 
                     filt.method, cutoff, trafo, 
                     log.n0, n.fold,
                     nsplit, ml.method, 
                     name.study = NULL, dir = NULL, 
                     max.show = 50, feat.type = "filtered", sd.add = 0.1, raw.auc = T, plots = T,
                     siamcat.holdout = NULL, use.holdout = F,
                     param.set = NULL,
                     cases = "PD", 
                     perform.selection = FALSE, 
                     param.selection = list(thres.fs = 100, method.fs = "AUC", direction="absolute")){
  require("SIAMCAT")
  meta<-meta[order(row.names(meta)),]
  df<-df[, order(names(df))]
  stopifnot(all(names(df) == rownames(meta)))
  sc.obj <- siamcat(feat=df, meta=meta, label=name.feat, case=cases)
  sc.obj <- filter.features(sc.obj,  
                            filter.method = filt.method,
                            cutoff = cutoff)
  sc.obj <- normalize.features(
    sc.obj,
    norm.method =trafo,
    norm.param = list(
      log.n0 = log.n0,
      n.p=2,
      norm.margin = 1,
      sd.min.q = sd.add),
    feature.type = feat.type
    )
  sc.obj <-  create.data.split(
    sc.obj,
    num.folds = n.fold, 
    num.resample = nsplit 
  )
  sc.obj <- train.model(
    sc.obj,
    method = ml.method,
    param.set = param.set,
    perform.fs = perform.selection,
    param.fs = param.selection
  )
  sc.obj <- make.predictions(sc.obj, siamcat.holdout = NULL)
  sc.obj <-  evaluate.predictions(sc.obj)
  if(plots == T) {
    model.evaluation.plot(sc.obj, fn.plot = paste0(dir, name.study, 'eval_plot_dnk.pdf'))
    
    model.interpretation.plot(
    sc.obj,
    fn.plot = paste0(dir, name.study, 'interpretation.pdf'),
    consens.thres = 0.5,
    limits = c(-3, 3),
    heatmap.type = 'zscore',
    max.show = max.show)
  }else {}
  perf<-c(sc.obj@eval_data$auroc, sc.obj@eval_data$auprc)
  names(perf)<-c("auroc", "auprc")
  if(raw.auc == T){
    return(list(perf, sc.obj))
  } else {return(sc.obj)}
}


#################################################################################
#
# Script to run SIAMCAT in parallel
# the number of cores indicate the fold the cross validation has to be repeated.
#
#################################################################################


siamcat.wf.parallel<-function(numCores = 5,  # these are the options for the parallel jobs
                              df, meta, name.feat, 
                              filt.method, cutoff, trafo, 
                              log.n0, n.fold = 10, # n.fold is the  cross validation n.fold
                              nsplit = 1, ml.method, 
                              name.study = NULL, dir = NULL, 
                              max.show = 50, feat.type = "filtered", sd.add = 0.1, raw.auc = T, 
                              plots = T, param.set = NULL, verbose = T, cases = "PD",  perform.selection = FALSE, 
                              param.selection = list(thres.fs = 100, method.fs = "AUC", direction="absolute"),
                              reps = reps) # verbose is to see what foreach does
  {
  require(foreach)
  require(doParallel)
  cl <- makeForkCluster(numCores)
  registerDoParallel(cl)
  set.seed(33) # for reproducibility so the list of seed is always the same even when seed is not set in the R session
  list.of.seeds<-sample(1:100000000, numCores*reps*1000, replace=FALSE)
  require(mlr3extralearners)
      if(reps > 1){
        f.f<-vector(mode = "list", length = reps)
        for(i in 1:reps){
          print(paste0("Running reps ", i, " of ", reps))
          print("")
          # print out the index and seeds used
          print(numCores)
          print(reps)
          nC<-1:numCores
          print("Index used to select seeds are: ")
          print(nC+(i^3))
          print("The seeds used are")
          print(list.of.seeds[nC+(i^3)])
          f.f[[i]]<-foreach (nC=1:numCores, .packages=c("SIAMCAT", "mlr3", "mlr3learners", "mlr3extralearners"), .verbose = verbose) %dopar% { 
            set.seed(list.of.seeds[nC+(i^3)])
            # useful to debug
            list.of.seeds[nC+(i^3)]
            nC+(i^3)

            siamcat.wf(df = df, meta = meta, name.feat = name.feat, 
                       filt.method= filt.method, cutoff=cutoff,
                       trafo=trafo, log.n0=log.n0,
                       ml.method=ml.method,
                       n.fold = n.fold,
                       nsplit =nsplit,
                       raw.auc = raw.auc,
                       plots = plots,
                       param.set = param.set, 
                       cases = cases, 
                       perform.selection = perform.selection,
                       param.selection = param.selection)
          }
          on.exit(stopCluster(cl))
        }
        return(f.f)
      } else {
        cores<-1:numCores
        print("Index used to select seeds are: ")
        print(cores)
        print("The seeds used are")
        print(list.of.seeds[cores])
        f<-foreach (nC=1:numCores, .packages=c("SIAMCAT", "mlr3", "mlr3learners", "mlr3extralearners"), .verbose = verbose) %dopar% { 
          set.seed(list.of.seeds[nC])
          siamcat.wf(df = df, meta = meta, name.feat = name.feat, 
                     filt.method= filt.method, cutoff=cutoff,
                     trafo=trafo, log.n0=log.n0,
                     ml.method=ml.method,
                     n.fold = n.fold,
                     nsplit =nsplit,
                     raw.auc = raw.auc,
                     plots = plots,
                     param.set = param.set, 
                     cases = cases, 
                     perform.selection = perform.selection,
                     param.selection = param.selection)
            }
            on.exit(stopCluster(cl))
        }
      return(f)
}



#################################################################################
#
# Script for calculating ROC on the output of the parallel ML calculation
# This script needs to be applied to a list of ML models (SIAMCAT obj)
#
#################################################################################

calculate.roc<-function(list.siamcat, 
                        meta,
                        auc.plot = T,
                        direction = "<"){
  tmp<-do.call(cbind.data.frame, lapply(list.siamcat, function(x) {
    x@pred_matrix
  }))
  tmp<-tmp[order(row.names(tmp)),]
  meta<-meta[order(rownames(meta)),]
  stopifnot(all(rownames(tmp) == rownames(meta)))
  tmp$average<-apply(tmp, 1, mean)
  meta$PD<-as.factor(meta$PD) 
  meta$PD<-relevel(meta$PD, ref = "HC")
   r<-pROC::roc(predictor = tmp$average, response = meta$PD, direction = direction) # direction indicate the ratios for comparison
 if(auc.plot == T){
    plot(r, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
         grid.col=c("green", "red"), max.auc.polygon=F,
         auc.polygon.col="lightblue", print.thres=F)  }
  return(r) 
}

#################################################################################
#
# Script for Building models on single dataset using parallel workflow
#
#################################################################################

siamcat.wf.lists<-function( vector.names, name.study.var, # vector that contains the names of the studies. It has to be exactly as is in the meta obj
  # below is the #cores for the paralel wf
  numCores = 5,  
  # below are the option for the siamcat.wf function
  df, meta, name.feat, 
  filt.method, cutoff, trafo, 
  log.n0, n.fold = 2,
  nsplit = 1, ml.method, 
  name.study = NULL, 
  dir = NULL, 
  max.show = 50, 
  feat.type = "filtered",
  sd.add = 0.1, raw.auc = F, 
  plots = F, param.set = NULL, 
  auc.plot = F, cases = "PD",
  save.all = T, # decide whether you want to save all objects or only the auroc
  reps = 1,  
  perform.selection = FALSE,
  param.selection = list(thres.fs = 100, 
                         method.fs = "AUC", 
                         direction="absolute")) 
  {
    print(param.selection)
  
    l_smct<-vector(mode = "list", length = length(vector.names)) # store for each study a list of models
    l_auroc<-vector(mode = "list", length = length(vector.names))

    names(l_smct)<-vector.names
    names(l_auroc)<-vector.names
    
    for(i in 1:length(vector.names)){
      print(vector.names[i])
      # select the relevant data
      tmp<-meta[ meta[,name.study.var]== vector.names[i],]
      print("")
      print("")
      print(paste0("Calculating models for ", vector.names[i]))
      tmp.feat<-df[,rownames(tmp)] 
      # remove the 0 taxa/samples
      tmp.feat <-tmp.feat[which(rowSums(tmp.feat) != 0), which(colSums(tmp.feat) != 0)]  
      if(length(which(rowSums(tmp.feat) == 0)) > 0){
        print(paste0("NOTE: we removed some rows: ", (length(which(rowSums(tmp.feat) == 0)))))
      }
      if(length(which(colSums(tmp.feat) == 0)) > 0){
        print(paste0("NOTE: we removed some Cols: ", (length(which(colSums(tmp.feat) == 0)))))
      }
      if(numCores == 0){
        # make a list to save the model NOT run in parallel
        l.not.parallel<-vector(mode = "list", length = length(vector.names))
        names(l.not.parallel)<-vector.names
        print("I will run the Siamcat_wf not in parallel. Be patient please!")
        require("mlr3")
        require("mlr3learners")
        require("mlr3extralearners")
        siam<-siamcat.wf(df = tmp.feat, meta = tmp, name.feat = name.feat, 
                                filt.method= filt.method, cutoff=cutoff,
                                trafo=trafo, log.n0=log.n0,
                                ml.method=ml.method,
                                n.fold = n.fold,
                                nsplit =nsplit,
                                raw.auc = raw.auc,
                                plots = plots,
                                param.set = param.set, 
                                cases = cases, 
                                perform.selection = perform.selection,
                                param.selection = param.selection)
        auc<-eval_data(siam)$auroc %>% as.double()
      } else {
        if(reps > 1){ # note: even if the numCors is 1 this will be run by the code below. 
          print(paste0("I am going to run the parallel workflow, to generate ", reps * numCores, " indipendent models, be patient!"))
          lsiam<-siamcat.wf.parallel(reps = reps, numCores = numCores, df = tmp.feat, meta = tmp, name.feat = name.feat,
                                         filt.method = filt.method,
                                         cutoff = cutoff, log.n0 = log.n0, nsplit = nsplit, n.fold = n.fold,
                                         trafo = trafo,  ml.method = ml.method,
                                         raw.auc = raw.auc, plots = plots,
                                         max.show = max.show, param.set = param.set, cases = cases,
                                         perform.selection = perform.selection,
                                         param.selection = param.selection)
          print("I am going to concatenate the lists now.")
          siam<-do.call(c, lsiam) # as the output of parallel is a list that needs to be concatenated
          print(paste0("The lenght of the parallel Siamcat list is ", length(siam)))
          auc<-calculate.roc(siam, meta = tmp, auc.plot = auc.plot)
          } else {
          print(paste0("I am going to run the parallel workflow, to generate ", reps * numCores, " indipendent models, be patient!"))
          siam<-siamcat.wf.parallel(reps = reps, numCores = numCores, df = tmp.feat, meta = tmp, name.feat = name.feat,  
                                          filt.method = filt.method, 
                                          cutoff = cutoff, log.n0 = log.n0, nsplit = nsplit, n.fold = n.fold, 
                                          trafo = trafo,  ml.method = ml.method, 
                                          raw.auc = raw.auc, plots = plots,
                                          max.show = max.show, param.set = param.set, cases = cases,
                                          perform.selection = perform.selection,
                                          param.selection = param.selection)
          print(paste0("The lenght of the parallel Siamcat list is ", length(siam)))
          auc<-calculate.roc(siam, meta = tmp, auc.plot = auc.plot)
          }
      }
      # save the list of models and AUC at the same if it was run in parallel
      if(save.all == T){
        print("Both SIAMCAT models and AUC will be saved!")
        l_smct[[i]]<-siam
        l_auroc[[i]]<-auc
      } else {
        print("Only AUC will be saved!")
        l_auroc[[i]]<-auc
        rm(siam)
        l_smct[[i]]<-NA
      }
     }
    
    # return the S4 object with results
    setClass("siamcat.wf.out", slots=list(siamcat="list", auc="list"))
    out<-new("siamcat.wf.out", siamcat = l_smct, auc = l_auroc)
    return(out)
}





#################################################################################
#
# Script for the LOSO validation
# This script does not work in parallel.
# 
#################################################################################

loso.validation<-function( vector.names, # this is a vector that contains the names of the studies. It has to be exactly as is in the meta object
                           # below are the option for the siamcat.wf function
                           df, meta,
                           name.feat,
                           filt.method, cutoff = 1e-07, trafo, 
                           log.n0, n.fold = 2, 
                           nsplit = 1, ml.method, 
                           name.study = NULL, dir = NULL,
                           max.show = 50, feat.type = "filtered", sd.add = 0.1, raw.auc = F, plots = F,
                           siamcat.holdout = NULL, param.set = NULL,
                           study.label = "Study",
                           case,
                           n.p=2,
                           norm.margin = 1, perform.selection = FALSE, 
                           param.selection = list(thres.fs = 100, method.fs = "AUC", direction="absolute")) 
{
  l_smct<-vector(mode = "list", length = length(vector.names)) # store for each study a list of models
  l_smct.test<-vector(mode = "list", length = length(vector.names)) # store for the test siamcat objects
  
  l_auroc<-vector(mode = "list", length = length(vector.names))
  names(l_smct)<-vector.names
  names(l_smct.test)<-vector.names
  names(l_auroc)<-vector.names
  
  df_loso_cv<-data.frame(study.out=vector.names, 
                         AUC = vector(length = length(vector.names)), 
                         n.samples =  vector(length = length(vector.names)),
                         n.pd = vector(length = length(vector.names)),
                         n.hc = vector(length = length(vector.names)))
  for(i in 1:length(vector.names)){
    print(paste0("Running model onf whole dataset leaving out ", vector.names[i]))
    # make training set
    tmp<-subset(meta, meta[, study.label] != vector.names[i])
    tmp.feat<-df[, rownames(tmp)] 
    
    # fill in info on sample size in df
    df_loso_cv[i,3]<-nrow(tmp)
    df_loso_cv[i,4]<-nrow(subset(tmp, PD == "PD"))
    df_loso_cv[i,5]<-nrow(subset(tmp, PD == "HC"))
    
    # create training siamcat object
    stopifnot(all(names(df) == rownames(meta)))
    siamcat.wf.parallel
    sc.obj <- siamcat(feat=tmp.feat, meta=tmp, label=name.feat, case=case)
    sc.obj <- filter.features(sc.obj,  
                              filter.method = filt.method,
                              cutoff = cutoff)
    sc.obj <- normalize.features(
      sc.obj,
      norm.method =trafo,
      norm.param = list(
        log.n0 = log.n0,
        n.p=n.p,
        norm.margin = norm.margin,
        sd.min.q = sd.add),
      feature.type = feat.type
    )
    sc.obj <-  create.data.split(
      sc.obj,
      num.folds = n.fold, 
      num.resample = nsplit
    )
    sc.obj <- train.model(
      sc.obj,
      method = ml.method,
      param.set = param.set,
      perform.fs = perform.selection,
      param.fs = param.selection
    )
    l_smct[[i]]<-sc.obj
    
    # make testing set
    tmp.test<-subset(meta,  meta[, study.label]  == vector.names[i])
    tmp.test.feat<-df[, rownames(tmp.test)] 
    
    # create Testing siamcat object
    sc.obj.test <- siamcat(feat=tmp.test.feat, meta=tmp.test,label= name.feat, case=case)
    # Running  model with holdout
    sc.obj.test <- make.predictions(sc.obj, siamcat.holdout = sc.obj.test)
    sc.obj.test <- evaluate.predictions(sc.obj.test)
    l_smct.test[[i]]<-sc.obj.test
    
    # extract auc
    df_loso_cv[i,2]<-eval_data(sc.obj.test)$auroc %>% as.double()
    
    
  }
  setClass("siamcat.wf.out", slots=list(siamcat.train="list", siamca.test = "list", auc="data.frame"))
  out<-new("siamcat.wf.out", siamcat.train = l_smct, siamca.test = l_smct.test, auc = df_loso_cv)
  return(out)
}
