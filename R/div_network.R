
#' Calculate a bipartite network for microbiome data
#'
#' @param ps phyloseq Object, contains OTU tables, tax table and map table, represented sequences,phylogenetic tree.
#' @param num Number of samples per group.
#' @param flour Is the arrangement of points discrete? TRUE or FEASE would be selected
#' @examples
#' data(ps)
#' ps_sub = filter_taxa(ps, function(x) sum(x ) > 20 , TRUE)
#' ps_sub = filter_taxa(ps_sub, function(x) sum(x ) < 30 , TRUE)
#' ps_sub
#'
#' result = div_network(ps_sub,num = 6)
#'
#' edge = result[[1]]
#' head(edge)
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export



div_network = function(ps,num = 6,group = "Group",flour = TRUE){
  mapping = as.data.frame(sample_data(ps))
  mapping = mapping[,group]
  colnames(mapping[,group]) <- "Group"
  sample_data(ps) = mapping
  ps_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps_rela
  vegan_otu <-  function(physeq){
    OTU <-  otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <-  t(OTU)
    }
    return(as(OTU,"matrix"))
  }
  vegan_tax <-  function(physeq){
    tax <-  tax_table(physeq)

    return(as(tax,"matrix"))
  }


  aa = vegan_otu(ps)
  otu_table = as.data.frame(t(aa))
  count = aa
  countA = count
  dim(count)
  sub_design <- as.data.frame(sample_data(ps))


  ##########这里的操作为提取三个分组
  pick_val_num <- num*2/3
  count[count > 0] <- 1###这个函数只能用于0,1 的数据，所以我这么转换

  count2 = as.data.frame(count)

  #数据分组
  iris.split <- split(count2,as.factor(sub_design$Group))
  #数据分组计算平均值
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  # 组合结果
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine)
  head(ven2)
  ven2[ven2 < pick_val_num]  = 0
  ven2[ven2 >=pick_val_num]  = 1
  ven2 = as.data.frame(ven2)


  ven3 = as.list(ven2)
  ven2 = as.data.frame(ven2)

  if (flour == TRUE) {
    ven2 = ven2[rowSums(ven2) == dim(ven2)[2] | rowSums(ven2) == 1,]

  }


  tax_table = as.data.frame(vegan_tax(ps))
  otu_table = as.data.frame(t(vegan_otu(ps)))
  dim(otu_table)


  # head(tax_table)
  otu_net = merge(ven2,tax_table,by = "row.names",all= F)

  row.names(otu_net) = otu_net$Row.names
  otu_net$Row.names = NULL
  head(otu_net)
  OTU = as.matrix(otu_table)
  norm = t(t(OTU)/colSums(OTU)) #* 100 # normalization to total 100
  norm1 = norm %>%
    t() %>% as.data.frame()

  #数据分组计算平均值
  iris.split <- split(norm1,as.factor(mapping$Group))

  iris.apply <- lapply(iris.split,function(x)colSums(x))
  # 组合结果
  norm2 <- do.call(rbind,iris.apply)%>%
    t()

  dim(norm2)
  colnames(norm2) = paste(colnames(norm2),"mean",sep = "")
  dim(otu_net)
  otu_net2 = merge(otu_net,norm2,by = "row.names",all= F)
  dim(otu_net2)
  head(otu_net2)
  colnames(otu_net2)[1] = "ID"
  # otu_net2$Label = paste("OTU_",1:length(otu_net2$ID),ep = "")
  ### 注释中缺乏样品标签，需要添加样品的注释信息
  sample_label = otu_net2[1:length(unique(mapping$Group)),]
  # sample_label = as.matrix(sample_label )


  sample_label$ID = unique(mapping$Group)
  # sample_label$Label = unique(mapping$SampleType)

  point_label = rbind(otu_net2,sample_label)

  #Kingdom+ Phylum + Class+ Order+ Family +Genus+Species+2Mmean+2MNPKmean +CKmean +MNPKmean + NPKmean ,value.var = "Abundance"
  # otu_net2$Label = NULL
  head(otu_net2)
  if (length(colnames(tax_table(ps))) == 6) {
    net_all = reshape2::melt(otu_net2, id=c("ID","Kingdom","Phylum" ,"Class", "Order","Family" ,"Genus",
                                            paste(unique(mapping$SampleType),"mean",sep = "")))
  }


  if (length(colnames(tax_table(ps))) == 7) {
    net_all = reshape2::melt(otu_net2, id=c("ID","Kingdom","Phylum" ,"Class", "Order","Family" ,"Genus","Species",
                                            paste(unique(mapping$Group),"mean",sep = "")))
  }
  if (length(colnames(tax_table(ps))) != 7|length(colnames(tax_table(ps))) != 6 ) {

    net_all = reshape2::melt(otu_net2, id=c("ID",rank_names(ps),
                                            paste(unique(mapping$Group),"mean",sep = "")))
  }

  net_filter<- dplyr::filter(net_all, value != 0 )


  net_fal = data.frame(
    source = net_filter$ID,
    target = net_filter$variable,
    connect = rep("pp",nrow(net_filter)),
    value = net_filter$value,
    Label =  net_filter$ID
  )
  head(net_fal)

  return(list(net_fal,point_label,ven2))

}

vegan_otu <-  function(physeq){
  OTU <-  otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <-  t(OTU)
  }
  return(as(OTU,"matrix"))
}

vegan_tax <-  function(physeq){
  tax <-  tax_table(physeq)

  return(as(tax,"matrix"))
}
