####################################################################
#
# Script to extract the coefficients from SIAMCAT models run in parallel
#
##################################################################

coefficient.extraction<-function(studies,
                                 list.models,
                                 model.parallel = T,
                                 string.to.add = F,
                                 select = "mean.weight", 
                                 return.all.coeff = F) {
  if(return.all.coeff == F) {
    if(model.parallel == T){
      list_coeff<-vector(mode = "list", length = length(studies))
      names(list_coeff)<-studies
      for(s in 1:length(studies)){
        # for this I run X different models in parallel. So there is a list of X siamcat objects and each has Y fold CV
        mod<-list.models@siamcat[[studies[s]]]
        l_coef<-lapply(mod, function(m) feature_weights(m))
        # the list above has X dfs with weights for the Y fold CV I need to make an average on this
        # calculate overall average
        df<-data.frame(taxa = rownames(l_coef[[1]]))
        for(i in 1: length(l_coef)){
          df[,paste0(select, "_", i)]<- l_coef[[i]][,select]
        }
        df$overall.average<-apply(df[,-1], 1, mean)
        if(string.to.add == F){
          df$study<-rep(studies[s], length = nrow(df))
        } else {df$study<-rep(paste0(studies[s], "_", string.to.add), length = nrow(df))}
        list_coeff[[s]]<-df
      }
    } else {
      list_coeff<-vector(mode = "list", length = length(studies))
      names(list_coeff)<-studies
      for(i in 1:length(studies)){
        mod<-list.models[[grep(paste0("^",studies[[i]]), names(list.models))[1]]]@siamcat.train 
        coef<-feature_weights(mod) # this is a result of X times Y fold CV
        if(string.to.add == F){
          coef$study<-rep(studies[i], length = nrow(coef))
        } else {coef$study<-rep(paste0(studies[i], "_", string.to.add), length = nrow(coef))}      
        list_coeff[[i]]<-coef
      }
    }
    return(list_coeff)
  } else {
    if(model.parallel == T){
      list_coeff<-vector(mode = "list", length = length(studies))
      names(list_coeff)<-studies
      
      list_coeff<-vector(mode = "list", length = length(studies))
      names(list_coeff)<-studies
      
      for(s in 1:length(studies)){
        # for this I run X different models in parallel. So there is a list of X siamcat objects and each has Y fold CV
        mod<-list.models@siamcat[[studies[s]]]
        # extract all models for all the Y CV in each X repetitions
        l.mod.unwr<-lapply(mod, function(x) x@model_list$models)
        
        #extract coeff - > they have all the same number of rows
        l.coeff<-lapply(l.mod.unwr, function(x) lapply(x, function(m) coef(getLearnerModel(m, more.unwrap = TRUE)) %>% 
                                                         as.matrix %>%
                                                         as.data.frame))
        # combine
        df.coeff<-do.call(cbind, do.call(cbind, l.coeff)) %>% as.data.frame()
        names(df.coeff)<-rep(paste0("W", 1:50))
        list_coeff[[s]]<-df.coeff
      }
      return(list_coeff)
    }
  }
}

  
  
  





