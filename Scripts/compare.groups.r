# Function to extract distances between samples
distance.between.groups <- function(dist.matrix, 
                                    meta, 
                                    variable.study = "Study_2", 
                                    same = FALSE) {
  require(dplyr)
  dist.m<-as.matrix(dist.matrix)
  stopifnot(isSymmetric(dist.m))
  studies <- unique(meta[,variable.study]) %>% 
    as.character()
  # make list to save items
  l.dist <- vector(mode = "list", length = length(studies))
  names(l.dist)<-studies
  for(id in studies) {
    print(id)
    study.same <- meta[which(meta[,variable.study] == id), "SampleID"]
    # double check
    if(same == FALSE){
      study.others<-meta$SampleID[-which(meta$SampleID %in% study.same)]
      stopifnot(all(!(study.same %in% study.others)))
    } else {study.others <- study.same
    stopifnot(all(study.same == study.others))
    }
    dist <- dist.m[study.same, study.others] 
    # double check
    if(same == FALSE){
      stopifnot(all(!(colnames(dist) %in% rownames(dist))))
      stopifnot(all(rownames(dist) == meta[which(meta[,variable.study] == id), "SampleID"]))
      
    } else {study.others <- study.same
    stopifnot(all(colnames(dist) == rownames(dist)))
    stopifnot(all(rownames(dist) == meta[which(meta[,variable.study] == id), "SampleID"]))
    
    }
    if(same){
      stopifnot(isSymmetric(dist))
      l.dist[[id]]<-dist[lower.tri(x = dist, diag = F)]
    } else { l.dist[[id]]<-as.vector(dist)}
    
  }
  return(l.dist)
}
