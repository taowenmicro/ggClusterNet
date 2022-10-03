

alpha = function(otu = NULL,
                 tax = NULL,
                 map = NULL,
                 ps = NULL,
                 group = "Group",
                 inde="Shannon",
                 sampling = TRUE,
                 Plot = TRUE){
  # path = "./alpha/"
  # dir.create(path)
  # library(phyloseq)
  # library(tidyverse)
  # library(vegan)
  # library(picante)      #picante 包加载时默认同时加载 vegan
  # library(agricolae)

  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)

  # map = as.data.frame(sample_data(ps))
  
  
  
  # ##按照最小序列数抽平
  # total = min(sample_sums(ps));total
  # # total = min(sample_sums(ps));total
  # standf = function(x,t = total)round(t*(x/sum(x)))
  # ps11 = transform_sample_counts(ps,standf)
  if (sampling == TRUE) {
    samplesize = min(phyloseq::sample_sums(ps))
    if (samplesize == 0) {
      print("0 number sequence of some samples")
      print("median number were used")
      ps11  = phyloseq::rarefy_even_depth(ps,sample.size = samplesize)
    } else{
      ps11  = phyloseq::rarefy_even_depth(ps,sample.size = samplesize)
    }
  } else if(sampling == FALSE){
    ps11 = ps
  }
  

  mapping = phyloseq::sample_data(ps11)
  ps11 = phyloseq::filter_taxa(ps11, function(x) sum(x ) >0 , TRUE); ps11
  head(mapping)
  colnames(mapping) = gsub(group,"AA", colnames(mapping))
  
  mapping$Group = mapping$AA
  mapping$Group = as.factor(mapping$Group)
  mapping$Group 

  count = as.data.frame(t(ggClusterNet::vegan_otu(ps11)))

  alpha=vegan::diversity(count, "shannon")
  
  x = t(count) ##转置，行为样本，列为OTU
  head(x)
  
  
  ##核心，计算两个指标
  Shannon = vegan::diversity(x)  ##默认为shannon
  Shannon
  Inv_Simpson <- vegan::diversity(x, index = "invsimpson")
  Inv_Simpson
  
  #计算OTU数量
  S <- vegan::specnumber(x);S  ##每个样本物种数。等价于S2 = rowSums(x>0)
  S2 = rowSums(x>0)
  
  
  #多样性指标：均匀度Pielou_evenness，Simpson_evenness
  Pielou_evenness <- Shannon/log(S)
  Simpson_evenness <- Inv_Simpson/S
  
  est <- vegan::estimateR(x)
  est <- vegan::estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  
  report = cbind(Shannon, Inv_Simpson, Pielou_evenness, Simpson_evenness,
                 Richness, Chao1,ACE) #不同列进行合并
  head(report)
  
  index = merge(mapping,report , by="row.names",all=F)

  return(index)
}






