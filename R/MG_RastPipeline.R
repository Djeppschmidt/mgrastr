#' download taxonomic data
#'
#' download taxonoic data for a single annotation run
#' @param sample mg-rast ID
#' @param auth mg-rast authentication code
#' @keywords mg-RAST taxonomy download
#' @export
#' @examples
#' downloadTaxa()

downloadTaxa<-function(sample, auth){
  require(httr)
  require(jsonlite)
  base<-"https://api-ui.mg-rast.org/metagenome/"
  end<-"?verbosity=stats&detail=taxonomy&auth="
  s.tax<-fromJSON(content(GET(paste(base,sample,end,auth,sep="")), "text"), flatten=T)
  s.tax
}

#' download taxonomic data
#'
#' this functio wraps and summarizes the download taxa so that many samples can be processed together
#' @param sample mg-rast ID, or list of IDs
#' @param auth mg-rast authentication code
#' @keywords mg-rast taxonomy download
#' @export
#' @examples
#' download.O()
download.O<-function(sample, auth){
  l.t<-lapply(sample, downloadTaxa, auth)
  l.t
}


#' extract domain
#'
#' extract domain level information from taxonomy mg-rast download
#' @param list list of samples
#' @keywords mg-RAST taxonomy download
#' @export
#' @examples
#' ag.Domain()
ag.Domain<-function(list){
  a<-list$domain
  a
}

#' Extract phylum
#'
#' description
#' @param list list of samples
#' @keywords mg-RAST taxonomy download
#' @export
#' @examples
#' ag.Phylum()
ag.Phylum<-function(list){
  a<-list$phylum
  a
}

#' Extract class
#'
#' description
#' @param list list of samples
#' @keywords mg-RAST taxonomy download
#' @export
#' @examples
#' ag.Class()
ag.Class<-function(list){
  a<-list$class
  a
}

#' Extract Order
#'
#' description
#' @param list list of samples
#' @keywords mg-RAST taxonomy download
#' @export
#' @examples
#' ag.Order()
ag.Order<-function(list){
  a<-list$order
  a
}

#' Extract family
#'
#' description
#' @param list list of samples
#' @keywords mg-RAST taxonomy download
#' @export
#' @examples
#' ag.family()
ag.family<-function(list){
  a<-list$family
  a
}

#' extract genus
#'
#' description
#' @param list list of samples
#' @keywords mg-RAST taxonomy download
#' @export
#' @examples
#' ag.Genus()
ag.Genus<-function(list){
  a<-list$genus
  a
}

#' extract species
#'
#' description
#' @param list list of samples
#' @keywords mg-RAST taxonomy download
#' @export
#' @examples
#' ag.Species()
ag.Species<-function(list){
  a<-list$species
  a
}

#' aggregate taxonomy table
#'
#' aggregate mg-rast taxonomy data across samples into table at specified taxonomy rank
#' @param list list of samples
#' @keywords mg-RAST taxonomy download
#' @export
#' @examples
#' ag.tax()
ag.tax<-function(l.taxa, f){
  require(reshape2)
  require(plyr)
  f<-ldply(l.taxa, f)# f is the function to be applied (taxanomic level)
  names(f)<-c("ID", "taxa", "values")
  t.f<-as.data.frame(dcast(f, taxa~ID))
  rownames(t.f)<-t.f$taxa
  t.f<-t.f[,-1]
  t.f[is.na(t.f)]<-0
  t.f[]<-lapply(t.f, as.numeric)
  t.f
}

#' condense
#'
#' combine annotations based on sample ID
#' @param df data frame of taxonomy or functional annotations
#' @param samkey sample key of unique mg-rast IDs in one column and sample IDs in another
#' @param sam sample metadata with rows being unique to each sample
#' @keywords mg-RAST taxonomy download
#' @export
#' @examples
#' condense()
condense<-function(df, samkey, sam){
  require(phyloseq)
  ps<-phyloseq(otu_table(df, taxa_are_rows = T), sample_data(samkey))
  ps2<-merge_samples(ps, "Sample_ID", fun = sum)
  sample_data(ps2)<-sam
  ps2
}


#' load
#'
#' submit function annotation query to mg-rast, get result url
#' @param x list of mg-rast accession numbers
#' @param auth authentication code for mg-rast
#' @param ont level of ontology. acceptable: Subsystems, COG, KO, NOG
#' @param level level of ontology. acceptable values: 1, 2, 3, 4.
#' @param E expected error value. default mg-rast is 5; (e^-5), larger numbers are more retrictive.
#' @param length minimun number of bp to match query string. mg-rast uses 15
#' @param id minimum percent match. mg-rast uses 60
#' @keywords mg-RAST function download
#' @export
#' @examples
#' load()
mg.load<-function(x, auth, ont, level, E, length, id){
  require(httr)
  require(jsonlite)
  end<-"&auth="
  s<-fromJSON(content(GET(paste("https://api.mg-rast.org//matrix/function?group_level=level", level, "&source=", ont, "&evalue=", E, "&identity=", id, "&length=", length, "&version=1&result_type=abundance&asynchronous=1&id=", x,"&auth=", auth, sep="")), "text"), flatten=T)
  u<-s$url
  names(u)<-names(x)
  u
}

#' load wrapper for many different samples
#'
#' load wrapper for many different samples
#' @param x list of mg-rast accession numbers
#' @param auth authentication code for mg-rast
#' @param ont type of ontology. acceptable: Subsystems, COG, KO, NOG
#' @param level level of ontology. acceptable values: 1, 2, 3, 4.
#' @param E expected error value. default mg-rast is 5; (e^-5), larger numbers are more retrictive.
#' @param length minimun number of bp to match query string. mg-rast uses 15
#' @param id minimum percent match. mg-rast uses 60
#' @param parallel logical; default is T
#' @keywords mg-RAST function download
#' @export
#' @examples
#' loadFunc()
loadFunc<-function(x, auth, ont, level, E, length, id, parallel){
  require(parallel)
  require(plyr)
  if(parallel==FALSE) {
    l<-lapply(x, mg.load, auth, ont, level, E, length, id)}
  else {
  l<-mclapply(x, mg.load, auth, ont, level, E, length, id)
  }
  l
}

#' download function annotation
#'
#' description
#' @param x list of mg-rast accession numbers
#' @param auth authentication code for mg-rast
#' @param ont level of ontology. acceptable: Subsystems, COG, KO, NOG
#' @keywords mg-RAST function download
#' @export
#' @examples
#' download.F()
download.F<-function(x, level, ont){
  require(httr)
  require(jsonlite)

  s.dl<-fromJSON(content(GET(x), "text"), flatten=T)

  if(ont=="Subsystems"){
    s.func<-matrix(data=NA, nrow=length(s.dl$data$rows$id),ncol=2)
    s.func[,2]<-s.dl$data$data
    if(level==3){
    s.func[,1]<-s.dl$data$rows$metadata.hierarchy.level3
    }
    else if(level==2){
      s.func[,1]<-s.dl$data$rows$metadata.hierarchy.level2
    }
    else if(level==1){
      s.func[,1]<-s.dl$data$rows$metadata.hierarchy.level1
    }
    else if(level==4){
     s.func[,1]<-s.dl$data$rows$id
    }
    else { return("error: level not specified, or out of bounds")}
  }

  else if (ont=="KO"){
    s.func<-matrix(data=NA, nrow=length(s.dl$data$rows$id),ncol=2)
    s.func[,2]<-s.dl$data$data

    if(level==3){
      s.func[,1]<-s.dl$data$rows$metadata.hierarchy.level3
    }
    else if(level==2){
     s.func[,1]<-s.dl$data$rows$metadata.hierarchy.level2
    }
    else if(level==1){
     s.func[,1]<-s.dl$data$rows$metadata.hierarchy.level1
    }
   else if(level==4){
     s.func[,1]<-s.dl$data$rows$id
    }
    else { return("error: level not specified, or out of bounds")}
  }

  else if (ont=="COG"){
    s.func<-matrix(data=NA, nrow=length(s.dl$data$rows$id),ncol=2)
    s.func[,2]<-s.dl$data$data

    if(level==3){
    s.func[,1]<-s.dl$data$rows$metadata.hierarchy.level3
    }
    else if(level==2){
      s.func[,1]<-s.dl$data$rows$metadata.hierarchy.level2
    }
    else if(level==1){
     s.func[,1]<-s.dl$data$rows$metadata.hierarchy.level1
    }
   else if(level==4){
    s.func[,1]<-s.dl$data$rows$id
    }
  else { return("error: level not specified, or out of bounds")}
  }

  else if (ont=="NOG"){
    s.func<-matrix(data=NA, nrow=length(s.dl$data$rows$id),ncol=2)
    s.func[,2]<-s.dl$data$data

    if(level==3){
      s.func[,1]<-s.dl$data$rows$metadata.hierarchy.level3
    }
    else if(level==2){
      s.func[,1]<-s.dl$data$rows$metadata.hierarchy.level2
    }
    else if(level==1){
      s.func[,1]<-s.dl$data$rows$metadata.hierarchy.level1
    }
    else if(level==4){
      s.func[,1]<-s.dl$data$rows$id
    }
  }

else {return("Error: ontology not specified, or out of bounds")}

  as.data.frame(s.func)
  s.func
}

#' download functions
#'
#' download functional annotations of many files, combine into a table
#' @param x list of mg-rast accession numbers
#' @param auth authentication code for mg-rast
#' @param ont level of ontology. acceptable: Subsystems, COG, KO, NOG
#' @keywords mg-RAST function download
#' @export
#' @examples
#' downloadFunc()
downloadFunc<-function(x, level, ont){
  require(reshape2)
  require(plyr)
  require(stats)
  #require(phyloseq)
  dl<-ldply(x, download.F, level, ont)
  names(dl)<-c("ID", "Function", "Values")
  dl$Values<-as.numeric(dl$Values)
  dl$Function<-as.character(dl$Function)
  dl<-na.omit(dl)
  #dl<-na.rm(dl)
  t.f<-as.data.frame(dcast(dl, Function~ID, fun.aggregate=sum, na.omit=T))
  rownames(t.f)<-t.f$Function
  t.f<-t.f[,-1]
  t.f[is.na(t.f)]<-0
  t.f[]<-lapply(t.f, as.numeric)
  t.f
} #works for Subsystems!!

