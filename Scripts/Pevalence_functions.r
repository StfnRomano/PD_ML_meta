# This is a function to infer "prevalence" of a taxa in a taxa-table
# There is also another function to filter taxa according with their prevalence in multiple data sets
# Prevalence here is defined as detection trashold and it does not take into account seq depth

prevalence.est<-function(df, taxa.are.rows = T){
  if(all(sapply(df, class) != "integer") & all(sapply(df, class) != "numeric")){
    stop("All data in the df need to be integers or numerics. Please reformat. =)")
  if(taxa.are.rows == F)
    stop("Taxa need to be rows. Please reformat, dude!")
  }
  require("dplyr")
  taxa<-rownames(df)
  df[df>0]<-1 # the df can be either counts or relative proportions
  samp<-ncol(df)
  prev<-(rowSums(df))/samp*100
  #names(prev)<-rownames(df)
  df.prev<-data.frame(taxa = taxa, prevalence = prev)
  return(df.prev)
}


prevalence.est.conditions<-function(df, taxa.are.rows = T, meta){
  if(all(sapply(df, class) != "integer") & all(sapply(df, class) != "numeric")){
    stop("All data in the df need to be integers or numerics. Please reformat. =)")
    if(taxa.are.rows == F)
      stop("Taxa need to be rows. Please reformat, dude!")
  }
  require("dplyr")
  taxa<-rownames(df)
  df[df>0]<-1 # the df can be either counts or relative proportions
  samp<-ncol(df)
  prev.tot<-(rowSums(df))/samp*100
  
  # PD
  pd<-subset(meta, PD =="PD")
  df.pd<-df[names(df) %in% pd$SampleID]
  samp<-ncol(df.pd)
  prev<-(rowSums(df.pd))/samp*100
  # HC
  hc<-subset(meta, PD =="HC")
  df.hc<-df[names(df) %in% hc$SampleID]
  samp<-ncol(df.hc)
  prev.hc<-(rowSums(df.hc))/samp*100
  #names(prev)<-rownames(df)
  df.prev<-data.frame(taxa = taxa, prevalence.hc = prev.hc, prevalence.pd = prev, prevalence = prev.tot)
  return(df.prev)
}



prevalence.matching<-function(list.df, n.studies = 2, prev = 5, taxa.are.rows = T){
  if(class(list.df) != "list"){
    stop("Input has to be a list of abundances etiher counts or perc, thanks dude!")
  }
  list.prev<-lapply(list.df, function(x) prevalence.est(x, taxa.are.rows = taxa.are.rows))
  require("dplyr")
  list.prev.df<-do.call(rbind.data.frame, list.prev) # create a df with all the data 
  tax<-unique(list.prev.df$taxa)
  tax.filt<-vector()
  for(i in 1:length(tax)){
    if(nrow(list.prev.df[list.prev.df$prevalence >= prev & list.prev.df$taxa == tax[i],]) >= n.studies){
      tax.filt<-c(tax.filt, tax[i])
    }
  }
  return(tax.filt)
}
