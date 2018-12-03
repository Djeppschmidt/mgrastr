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
  if(level==4){
  s<-fromJSON(content(GET(paste("https://api.mg-rast.org//matrix/function?group_level=function&source=", ont, "&evalue=", E, "&identity=", id, "&length=", length, "&version=1&result_type=abundance&asynchronous=1&id=", x,"&auth=", auth, sep="")), "text"), flatten=T)}
  else{s<-fromJSON(content(GET(paste("https://api.mg-rast.org//matrix/function?group_level=level", level, "&source=", ont, "&evalue=", E, "&identity=", id, "&length=", length, "&version=1&result_type=abundance&asynchronous=1&id=", x,"&auth=", auth, sep="")), "text"), flatten=T)}
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
  s.func # restructure to return a phyloseq w/ tax table extracted
}


#####
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
  t.f<-as.data.frame(dcast(dl, Function~ID, fun.aggregate=sum, na.omit=T))
  rownames(t.f)<-t.f$Function
  t.f<-t.f[,-1]
  t.f[is.na(t.f)]<-0
  t.f[]<-lapply(t.f, as.numeric)
  t.f
} #works for Subsystems!!

#' combine tables and condense
#'
#' download functional annotations of many files, combine into a table
#' @param x list of mg-rast accession numbers
#' @param auth authentication code for mg-rast
#' @param ont level of ontology. acceptable: Subsystems, COG, KO, NOG
#' @keywords mg-RAST function download
#' @export
#' @examples
#' DFTax()
DFTax<-function(l){
  l.t<-lapply(l, extractT)
  t<-ldply(l.t, cbind)
  t<-t[!duplicated(t$ID),]
  #t<-t[,2:6]
  rownames(t)<-t$ID
  t
}


#' parse single ontology into heirarchical table
#'
#' download functional annotations of many files, combine into a table
#' @param x list of mg-rast accession numbers
#' @param auth authentication code for mg-rast
#' @param ont level of ontology. acceptable: Subsystems, COG, KO, NOG
#' @keywords mg-RAST function download
#' @export
#' @examples
#'extractT()
extractT<-function(url){
  require(jsonlite)
  require(httr)
  u<-fromJSON(content(GET(url), "text"), flatten=T)
  t<-data.frame("L1"=u$data$rows$metadata.hierarchy.level1,"L2"=u$data$rows$metadata.hierarchy.level2,"L3"=u$data$rows$metadata.hierarchy.level3,"L4"=u$data$rows$metadata.hierarchy.level4, "ID"=u$data$rows$id)
  t
}

#' parse single ontology into heirarchical table
#'
#' download functional annotations of many files, combine into a table
#' @param x list of mg-rast accession numbers
#' @param auth authentication code for mg-rast
#' @param ont level of ontology. acceptable: Subsystems, COG, KO, NOG
#' @keywords mg-RAST function download
#' @export
#' @examples
#'dendromap()
dendromap<-function(){

}


#' construct table of network index values
#'
#' download functional annotations of many files, combine into a table
#' @param x list of phyloseq objects
#' @param groups1 list of groups to divide phylsoeq object by
#' @param groups2 groups for stats model
#'
#' @export
#' @examples
#'netstat()
netStat<-function(list, groups2, name){
  require(plyr)
  require(dplyr)
  #pdf(paste("~/Desktop/PhD/Metagenome/Networks/Disease.pdf", sep=""))
 # y$membership<-
  l_ply(list, plotNtwk, name)
  #while (!is.null(dev.list()))  dev.off()
  y$metrics<-ldply(list, ConnStat)
  y$metrics<-data.frame(y$metrics, groups2) # sanity check
  y$stats<-c("Closeness"=summary(aov(Closeness~Codes, data=y$metrics)),
             "Degree"=summary(aov(Closeness~Codes, data=y$metrics)),
             "Modularity"=summary(aov(Closeness~Codes, data=y$metrics)))
  y
}

#' set filter for function abundance
#'
#' download functional annotations of many files, combine into a table
#' @param b functions
#' @export
#' @examples
#'netstat()
filt<-function(b){sum(b>2)>(0.8*length(b))}

#' construct table of connectivity values
#'
#' download functional annotations of many files, combine into a table
#' @param x phyloseq object
#' @param cat categorical variable to subset phyloseq object
#'
#' @keywords igraph metric summaries
#' @export
#' @examples
#' netstat()
ConnStat<-function(list){
  require(phyloseq)
  require(igraph)
  #groups1<-groups1
  #a<-subset_samples(x, samplecodes==groups1)
  a<-filter_taxa(list, filt, TRUE)
  a<-transform_sample_counts(a, transform)
  a<-otu_table(a)
  a<-cor(a)
  a[abs(a)<0.7]<-0
  a[abs(a)>0.7]<-1
  net<-graph_from_incidence_matrix(a)
  cfg<- cluster_fast_greedy(as.undirected(net))
  out<-matrix(1:3,1)
  colnames(out)<-c("Closeness", "Degree", "Modularity")
  out[1,1]<-mean(closeness(net))
  out[1,2]<-mean(degree(net))
  out[1,3]<-modularity(net, membership(cfg))
  out
}

#' construct table of connectivity values
#'
#' download functional annotations of many files, combine into a table
#' @param x phyloseq object
#' @param cat categorical variable to subset phyloseq object
#'
#' @keywords igraph metric summaries
#' @export
#' @examples
#'plotNtwk()
plotNtwk<-function(list, name){
  require(phyloseq)
  require(igraph)
  #a<-subset_samples(x, samplecodes==groups1)
  a<-filter_taxa(list, filt, TRUE)
  a<-transform_sample_counts(a, transform)
  a<-otu_table(a)
  a<-cor(a)
  a[a<0.7]<-0
  a[a>0.7]<-1
  net<-graph_from_incidence_matrix(a)
  cfg<- cluster_fast_greedy(as.undirected(net))
  #ceb <- cluster_edge_betweenness(net)
  #out$plot1<-dendPlot(ceb, mode="hclust")
  #png(paste("~/Desktop/PhD/Metagenome/Networks/", x, ".png", sep=""))
  plot(cfg, as.undirected(net), layout=layout_nicely(net), vertex.label=NA, main=name)
  #dev.off()
  #out<-membership(cfg)
  out}

#' transform count values
#' @param x taxa
#' @keywords igraph metric summaries
#' @export
#' @examples
#'transform()
transform<-function(x){x/sum(x)}
