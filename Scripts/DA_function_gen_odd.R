# This function runs genOdds for case control studies
# It requires:
# * a list of df containing the data. df in the list have to have row = sample names; col = feats
# * the meta data matching the samples in the list above. Need a variable called Cohort and one called Sample_ID

calcualte_da_genOdds<-function( list.df,
                           meta,
                           treat = "PD",
                           assume_no_effect = FALSE){
  require(dplyr)
  require(genodds)
  # create list with results
  l.df<-list.df
  # loop through all the df in the list
  for(d in 1:length(list.df)) {
    cohort<-names(list.df)[[d]]
    print(paste0("Testing for ", cohort))
    tmp.spec<-list.df[[cohort, drop = F]]

    df<-data.frame( bin = vector(length = ncol(tmp.spec)),
                         estimate = vector(length = ncol(tmp.spec)),
                         std.error = vector(length = ncol(tmp.spec)),
                         p.value = vector(length = ncol(tmp.spec)))
   
    meta.tmp<-meta[meta$Cohort == cohort, ]
    # order meta and spec so they are in the same order and they can be merged
    meta.tmp<-meta.tmp[order(meta.tmp$Sample_ID),] 
    tmp.spec<-tmp.spec[order(rownames(tmp.spec)), ,drop =F]
  
    for(s in 1:ncol(tmp.spec)) {
      spec<-names(tmp.spec)[s] #  take naem of sepcies  
      print(paste0("This is variable number ", s, ", called: ", spec))
      
      if(all(meta.tmp$Sample_ID != rownames(tmp.spec))) {
        stop("Metadata and Spec abund table are not in the same order")
      } else {
        spec.ab<-tmp.spec[,s]
        tmp<-cbind.data.frame(spec.ab, meta.tmp) # build df with species and metadata
        x<-try(lo<-genodds(tmp$spec.ab,
                           tmp[,treat], 
                           assume_no_effect = assume_no_effect)) # in case of error
        if(class(x)[1] == "try-error"){
          df[s,"bin"]<-spec
          df[s,"estimate"]<-NA
          df[s,"std.error"]<-NA
          df[s,"p.value"]<-NA
        } else {
          df[s,"bin"]<-spec
          df[s,"estimate"]<-lo$pooled_lnodds
          df[s,"std.error"]<-lo$pooled_SElnodds
          df[s,"p.value"]<-lo$pooled_p
        }
    }
  }
    l.df[[d]]<-df
  }
  # now I need to bind the data frames together
  df.merged<-do.call(rbind, l.df)
  # add cohort to df
  df.merged$Cohort<-gsub("[.].*", "", rownames(df.merged))
  
  return(df.merged)
}

