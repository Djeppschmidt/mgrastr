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
  l_ply(list, plotNtwk, name)
  y<-NULL
  y$metrics<-ldply(list, ConnStat)
  y$metrics<-data.frame(y$metrics, groups2) # sanity check
  y$stats<-c("Closeness"=summary(aov(Cohesion~Codes, data=y$metrics)),
             "Degree"=summary(aov(Degree~Codes, data=y$metrics)),
             "Modularity"=summary(aov(Modularity~Codes, data=y$metrics)),
             "Nodes"=summary(aov(Nodes~Codes, data=y$metrics)),
             "Edges"=summary(aov(Edges~Codes, data=y$metrics)))
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
#' constructs network based explicitly on number of expected edges.
#' Challenge is in defining appropriate level of edges for interpretable topology.
#' Most nodes have at least one highly significant correlation to another; thus with very few edges, we get all the nodes, and as we decrease significance threshold, the network fills out.
#' @param x phyloseq object
#' @param cat categorical variable to subset phyloseq object
#'
#' @keywords igraph metric summaries
#' @export
#' @examples
#' netstat()
ConnStat2<-function(list, num){
  require(phyloseq)
  require(igraph)
  a<-filter_taxa(list, filt, TRUE)
  a<-transform_sample_counts(a, transform)
  o<-otu_table(a)
  c<-cor(o)
  i=1
  c[c<1]<-0
  n<-graph_from_incidence_matrix(c)
  while(ecount(n)<num){
    t<-otu_table(a)
    t<-cor(t)
    t[t<i]<-0
    t[t>i]<-1
    n<-graph_from_incidence_matrix(t)
    i=i-0.001
  }
  cfg<-cluster_fast_greedy(as.undirected(n))
  out<-matrix(1:3,1)
  colnames(out)<-c("Closeness", "Degree", "Modularity")
  out[1,1]<-mean(closeness(n))
  out[1,2]<-mean(degree(n))
  out[1,3]<-modularity(n, membership(cfg))
  out
}

#' plot networks
#'
#' make network explicitly defined number of edges
#' @param x phyloseq object
#' @param cat categorical variable to subset phyloseq object
#'
#' @keywords igraph metric summaries
#' @export
#' @examples
#'plotNtwk()
plotNtwk2<-function(list, name, num){
  require(phyloseq)
  require(igraph)
  a<-filter_taxa(list, filt, TRUE)
  a<-transform_sample_counts(a, transform)
  o<-otu_table(a)
  c<-cor(o)
  i=1
  c[c<1]<-0
  n<-graph_from_incidence_matrix(c)
  while(ecount(n)<num){
    t<-otu_table(a)
    t<-cor(t)
    t[t<i]<-0
    t[t>i]<-1
    n<-graph_from_incidence_matrix(t)
    i=i-0.001
  }
  cfg<-cluster_fast_greedy(as.undirected(n))
  plot(cfg, as.undirected(n), layout=layout_nicely(n), vertex.label=NA, main=name, vertex.size=10)
  #dev.off()
  #out<-membership(cfg)
  out}

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
  #filter for a priori relativized datasets (wrong approach!):
  #a<-filter_taxa(list, function(x) var(x) > 1e-5, TRUE)
  a<-transform_sample_counts(a, transform)
  o<-otu_table(a)
  c<-cor(o)
  c[abs(c)<0.7]<-0
  c[abs(c)>0.7]<-1
  n<-graph_from_incidence_matrix(c)
  cfg<-cluster_fast_greedy(as.undirected(n))
  out<-matrix(nrow=1, ncol=5)
  colnames(out)<-c("Cohesion", "Degree", "Modularity", "Nodes", "Edges")
  out[1,1]<-cohesion(n)
  out[1,2]<-mean(degree(n))
  out[1,3]<-modularity(n, membership(cfg))
  out[1,4]<-vcount(n)
  out[1,5]<-ecount(n)
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
  a<-filter_taxa(list, filt, TRUE)
  #filter for a priori relativized datasets (wrong approach!):
  #a<-filter_taxa(list, function(x) var(x) > 1e-5, TRUE)
  #a<-transform_sample_counts(a, transform)
  o<-otu_table(a)
  c<-cor(o)
  c[abs(c)<0.7]<-0
  c[abs(c)>0.7]<-1
  n<-graph_from_incidence_matrix(c)
  cfg<-cluster_fast_greedy(as.undirected(n))
  plot(cfg, as.undirected(n), layout=layout_nicely(n), vertex.label=NA, main=name, vertex.size=10)
  #dev.off()
  out<-membership(cfg)
  out}




#' transform count values
#' @param x taxa
#' @keywords igraph metric summaries
#' @export
#' @examples
#'transform()
transform<-function(x){x/sum(x)}

#' transform count values
#' @param x taxa
#' @keywords igraph metric summaries
#' @export
#' @examples
#'MeanNetval()
MeanNetval<-function(x){
  Means<-matrix(nrow=6, ncol=4)
  Means[1,1]<-mean(x$metrics$Modularity[x$metrics$Codes=="Reference"])
  Means[1,2]<-mean(x$metrics$Modularity[x$metrics$Codes=="Remnant"])
  Means[1,3]<-mean(x$metrics$Modularity[x$metrics$Codes=="Turf"])
  Means[1,4]<-mean(x$metrics$Modularity[x$metrics$Codes=="Ruderal"])


  Means[2,1]<-mean(x$metrics$Degree[x$metrics$Codes=="Reference"])
  Means[2,2]<-mean(x$metrics$Degree[x$metrics$Codes=="Remnant"])
  Means[2,3]<-mean(x$metrics$Degree[x$metrics$Codes=="Turf"])
  Means[2,4]<-mean(x$metrics$Degree[x$metrics$Codes=="Ruderal"])

  Means[3,1]<-mean(x$metrics$Closeness[x$metrics$Codes=="Reference"])
  Means[3,2]<-mean(x$metrics$Closeness[x$metrics$Codes=="Remnant"])
  Means[3,3]<-mean(x$metrics$Closeness[x$metrics$Codes=="Turf"])
  Means[3,4]<-mean(x$metrics$Closeness[x$metrics$Codes=="Ruderal"])


  Means[4,1]<-var(x$metrics$Modularity[x$metrics$Codes=="Reference"])
  Means[4,2]<-var(x$metrics$Modularity[x$metrics$Codes=="Remnant"])
  Means[4,3]<-var(x$metrics$Modularity[x$metrics$Codes=="Turf"])
  Means[4,4]<-var(x$metrics$Modularity[x$metrics$Codes=="Ruderal"])


  Means[5,1]<-var(x$metrics$Degree[x$metrics$Codes=="Reference"])
  Means[5,2]<-var(x$metrics$Degree[x$metrics$Codes=="Remnant"])
  Means[5,3]<-var(x$metrics$Degree[x$metrics$Codes=="Turf"])
  Means[5,4]<-var(x$metrics$Degree[x$metrics$Codes=="Ruderal"])

  Means[6,1]<-var(x$metrics$Closeness[x$metrics$Codes=="Reference"])
  Means[6,2]<-var(x$metrics$Closeness[x$metrics$Codes=="Remnant"])
  Means[6,3]<-var(x$metrics$Closeness[x$metrics$Codes=="Turf"])
  Means[6,4]<-var(x$metrics$Closeness[x$metrics$Codes=="Ruderal"])
  rownames(Means)<-c("ModularityM", "DegreeM", "ClosenessM","ModularityV", "DegreeV", "ClosenessV")
  colnames(Means)<-c("Reference", "Remnant", "Turf", "Ruderal")
  as.data.frame(Means)
}

#' Linear model
#' @param x taxon
#' @param y covariate data
#' @param data data
#' @keywords
#' @export
#' @examples
#'linmod()
linmod<-function(x, y, data){
  nlme(x~y, random=~1|City/Codes, data=data)
}

#' function-correlation matrix
#' @param x taxon
#' @param y functional abundance
#' @param data data
#' @keywords
#' @export
#' @examples
#'linmod()
FT.cor<-function(x, y){
  F.L4<-x
  F.L3<-tax_glom(x, taxrank = "L3")
  F.L1<-tax_glom(x, taxrank = "L1")

  T.a<-y
  T.g<-tax_glom(y, taxrank = "Genus")
  T.f<-tax_glom(y, taxrank = "Family")
  T.o<-tax_glom(y, taxrank = "Order")
  l<-c(as.data.frame(as.matrix(otu_table(T.a))), as.data.frame(as.matrix(otu_table(T.g))), as.data.frame(as.matrix(otu_table(T.f))), as.data.frame(as.matrix(otu_table(T.o))))

  list<-c("L1"=llply(l, cor, as.data.frame(as.matrix(otu_table(F.L1))), method="pearson"), "L3"=llply(l, cor, as.data.frame(as.matrix(otu_table(F.L3))), method="pearson"), "L4"=llply(l, cor, as.data.frame(as.matrix(otu_table(F.L4))), method="pearson"))
  list
}

#' function-correlation matrix
#' @param x taxon
#' @param y functional abundance
#' @param data data
#' @keywords
#' @export
#' @examples
#'linmod()

