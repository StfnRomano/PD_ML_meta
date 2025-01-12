
cross.study.validation<-function(list.models,
                                 microbiome,
                                 meta,
                                 names.studies,
                                 label="PD", 
                                 case="PD",
                                 auc.plot = F,
                                 name.study.var = "Study",
                                 study.in.models = T,
                                 meta.test.hold = NULL,
                                 microbiome.test.hold = NULL,
                                 name.study.var.hold = NULL,
                                 label.hold = NULL,
                                 case.hold = NULL,
                                 direction = "<",
# is.model.list ask whether the models for each studies are organized in a list of models. This is to accommodate the data structure resulting from the parallel SIAMCAT runs.
# The results of these runs is a siamcat object, which is a list of length = #Studies. Each study than contains the X number of models as a result of reps x ncores
                                 
                                 is.model.list = FALSE){ 
  if(is.model.list == TRUE){
    
  # take the name of the studies
    if(study.in.models == T){
      std<-names(list.models@siamcat)
      stopifnot(all(std == names.studies))
      # get the list of test sets
      l_csv<-vector(mode = "list", length = length(std))
      names(l_csv)<-std
  
      for(i in 1:length(std)){
        print(paste("The model to test is", std[i]))
        # get the name of the first model
        train<-std[i]
      # define the test data to use for this model
        test.data<-names.studies[-1*which(names.studies == train)]
      # select the list of models to test
        l_train<-list.models@siamcat[std[i]][[1]]
        l_all_tests<-vector(mode = "list", length = length(test.data))
      # run cross study validation
        for(tst in 1:length(test.data)){
          print(paste("Testing the model on", test.data[tst]))
          test.study<-test.data[tst]
        # select data that need to be tested
          meta.test<-meta[meta[, name.study.var] == test.study,]
          microbiome.test<-microbiome[, meta.test$SampleID]
        # create a SIAMCAT obj with that
          sc.obj.test <- siamcat(feat=microbiome.test, meta=meta.test,label=label, case=case)
        # make prediction with each model in the list of training data
          l_prediction <- lapply(l_train, function(x) {
            make.predictions(x, siamcat.holdout = sc.obj.test)
          })
          roc<-calculate.roc(l_prediction, meta.test, auc.plot = auc.plot, direction=direction)
        # prepare a vector containing the details of the tests
          out<-c(train, test.study, as.double(roc$auc))
          names(out)<-c("training.data", "testing.data", "AUC")
        # this is the list that contains the vector, roc obs, and list of 5 test sets
          l_test<-list(out, roc, l_prediction)
          names(l_test)<-c("summary", "roc", "testing.siamcat")
          print("Appending the resulting lists to final list")
          l_all_tests[[tst]]<-l_test
          print("Rename list")
          names(l_all_tests)[[tst]]<-paste0(std[i], "_vs_", test.study)
        }
        l_csv[[i]]<-l_all_tests
      }
      return(l_csv)
    } else {
      std<-names(list.models@siamcat) 
    # get the list of test sets
      l_csv<-vector(mode = "list", length = length(std))
      names(l_csv)<-std
    
      for(i in 1:length(std)){
        print(paste("The model to test is", std[i]))
        # get the name of the first model
        train<-std[i]
        # define the test data to use for this model
        test.data<-names.studies
        # select the list of models to test
        l_train<-list.models@siamcat[std[i]][[1]]
        l_all_tests<-vector(mode = "list", length = length(test.data))
        # run cross study validation
        for(tst in 1:length(test.data)){
          print(paste("Testing the model on", test.data[tst]))
          test.study<-test.data[tst]
          # select data that need to be tested
          meta.test<-meta.test.hold[meta.test.hold[, name.study.var.hold] == test.study,] # these are other meta for the testing data
          microbiome.test<-microbiome.test.hold[, meta.test$SampleID] # these are other microbiome for the testing data
          # create a SIAMCAT obj with that
          sc.obj.test <- siamcat(feat=microbiome.test, meta=meta.test,label=label.hold, case=case.hold)
          # make prediction with each model in the list of training data
          l_prediction <- lapply(l_train, function(x) {
            make.predictions(x, siamcat.holdout = sc.obj.test)
          })
          
          roc<-calculate.roc(l_prediction, meta.test, auc.plot = auc.plot, direction=direction) 
          # prepare a vector containing the details of the tests
          out<-c(train, test.study, as.double(roc$auc))
          names(out)<-c("training.data", "testing.data", "AUC")
          # this is the list that contains the vector, roc obs, and list of 5 test sets
          l_test<-list(out, roc, l_prediction)
          names(l_test)<-c("summary", "roc", "testing.siamcat")
          print("Appending the results list to final list")
          l_all_tests[[tst]]<-l_test
          print("Rename list")
          names(l_all_tests)[[tst]]<-paste0(std[i], "_vs_", test.study)
        }
        l_csv[[i]]<-l_all_tests
      }
      return(l_csv)
    }
  } else {
        
      std<-names(list.models)
      # get the list of test sets
      l_csv<-vector(mode = "list", length = length(std))
      names(l_csv)<-std
      
      for(i in 1:length(std)){
        print(paste("The model to test is", std[i]))
        # get the name of the first model
        train<-std[i]
        # define the test data to use for this model
        test.data<-names.studies
        # select the list of models to test
        l_train<-list.models[[i]]@siamcat.train
        l_all_tests<-vector(mode = "list", length = length(test.data))
        # run cross study validation
        for(tst in 1:length(test.data)){
          print(paste("Testing the model on", test.data[tst]))
          test.study<-test.data[tst]
          # select data that need to be tested
          meta.test<-meta.test.hold[meta.test.hold[, name.study.var.hold] == test.study,] # these are other meta for the testing data
          microbiome.test<-microbiome.test.hold[, meta.test$SampleID] # these are other microbiome for the testing data
          # create a SIAMCAT obj with that
          sc.obj.test <- siamcat(feat=microbiome.test, meta=meta.test,label=label.hold, case=case.hold)
          # make prediction with each model in the list of training data
          prediction <- make.predictions(l_train, siamcat.holdout = sc.obj.test)
        
          roc<-evaluate.predictions(prediction) 
          # prepare a vector containing the details of the tests
          out<-c(train, test.study, as.double(roc@eval_data$auroc))
          names(out)<-c("training.data", "testing.data", "AUC")
          # this is the list that contains the vector, roc obs, and list of 5 test sets
          l_test<-list(out, prediction)
          names(l_test)<-c("summary", "testing.siamcat")
          print("Appending the results list to final list")
          l_all_tests[[tst]]<-l_test
          print("Rename list")
          names(l_all_tests)[[tst]]<-paste0(std[i], "_vs_", test.study)
        }
        l_csv[[i]]<-l_all_tests
      }
      return(l_csv)
      }
    
}

###########################################################################
#
#
# This scripts extract the threshold corresponding to the 10 FPR form the training data
# Then it select all prediction in the testing above that threshold
# Divide this number by the number of samples in the testing data
#
###################################################################


cross.disease.pred<-function(list.of.csv.models,
                             list.train.models) {
   require(SIAMCAT)
  require(pROC)
  # do some check if they are not in the same order
    train.models<-names(list.of.csv.models)
    test<-names(list.train.models@siamcat)
    stopifnot(all(train.models==test))
    # create a dataframe in which to store the results
    empty<-c()
    for(i in 1:length(train.models)){
      training.data<-train.models[i]
      # extract the roc for each model
      t.auc<-list.train.models@auc[[i]]
      size.train<-list.train.models@siamcat[[i]][[1]]@label$label %>% length()
      train.auc<-t.auc$auc %>% as.double()
     ##################################
     # Extract the probability limit
     #######################################
      sp<-t.auc$specificities # here we extract the specificity from the auc predictions
      m<-which.min(abs(sp - 0.9))
      fpr.value<-abs(1-sp)[m] # get the index that we are actually selecting
      # extract the threshold that correspond to 10% FPR
      tr<-t.auc$thresholds[m]
      threshold.train<-tr
      # now use this threshold to see how many samples were predicted above this in the testing data set of the other disease
      mod<-list.of.csv.models[[i]]
      for(t in 1:length(mod)){
        train_vs_test<-names(mod)[t]
        testing.data<-gsub(".*_vs_", "", names(mod)[t])
        pro.test<-lapply(list.of.csv.models[[i]][[t]]$testing.siamcat, function(x) pred_matrix(x)) 
        # the above gives me a prediction probability for all samples for all X repetition of the Y models in the CV.
        # So in total I have X*Y probability per sample
        # make a df and calculate the average
        pro.test.df<-do.call(cbind, pro.test)
        prot.average<-apply(pro.test.df, 1, mean)
        x<-prot.average[prot.average >= tr] # select all probabilities >= to the 10% FPR threshold
        # calculate prediction rates (ration of samples > threshold and total sample in testing data)
        pr<-(length(x)/length(prot.average))*100 # prob vector contains all values for all samples in the testing set. So has the same length of testing set.
        
        tmp<-c(training.data, 
               size.train,
               testing.data,
               train.auc,
               train_vs_test,
               fpr.value,
               m, 
               threshold.train,
               pr)
        empty<-rbind.data.frame(empty, tmp)
      }
    }
    names(empty)<-c("training.data",
                    "size.train",
                    "testing.data" ,
                    "train.auc",
                    "train_vs_test",
                    "fpr.value",
                    "index.fpr",
                    "threshold.train",
                    "perc.pred")
    return(empty)
}

