#
# j = "Phylum"
# j = "Class"
# j = "Order"
# j =  "Family"
# j = "Genus"
# 
# otu =NULL
# tax = NULL
# map = NULL
# ps = ps1
# 
# group = "Group"
# axis_ord = NULL
# label = FALSE
# sd = FALSE
# Top = 10
# 
# tran = TRUE

# library(phyloseq)
# library(tidyverse)
# library(vegan)
# library(reshape2)
# library("plyr")
library(ggalluvial)
# detach("package:ggalluvial")
# detach("package:plyr")
# detach("package:reshape2")
# library(ggplot2)



barMainplot = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group",
                       j = "Phylum",axis_ord = NULL,label = TRUE ,sd = FALSE,Top = 10,tran = TRUE){
  
  
  if (is.null(axis_ord)) {
    axis_order = NA
  }else{
    axis_order = strsplit(basename(axis_ord), "-")[[1]]
  }
  
  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  
  # psdata = ps %>%
  #   tax_glom(taxrank = j)
  psdata <- ggClusterNet::tax_glom_wt(ps = ps,ranks = j)
  
  # transform to relative abundance
  if (tran == TRUE) {
    psdata = psdata%>%
      phyloseq::transform_sample_counts(function(x) {x/sum(x,na.rm = TRUE)} )
  }
  
  
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  phyloseq::tax_table(psdata)= tax
  
  Taxonomies <- psdata %>% # Transform to rel. abundance
    phyloseq::psmelt()
  

  Taxonomies$Abundance = Taxonomies$Abundance * 100
  # Taxonomies$Abundance = Taxonomies$Abundance /rep
  colnames(Taxonomies) <- gsub(j,"aa",colnames(Taxonomies))
  data = c()
  i = 2
  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]

    b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]

    c <- Taxonomies %>% 
      dplyr::filter(Group == a)
    c$Abundance <- c$Abundance/b
    data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)

    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }

  }
  # head(table)
  # sum(table$Abundance)
  # Taxonomies$Abundance = NULL
  # Taxonomies <- Taxonomies %>% inner_join(table,by = "Sample")
  # head(Taxonomies)
  Taxonomies = table


  
  #按照分组求均值

  by_cyl <- dplyr::group_by(Taxonomies, aa,Group)

  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance,na.rm = TRUE), sd(Abundance,na.rm = TRUE))

  
  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance,na.rm = TRUE))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))
  
  ##使用这个属的因子对下面数据进行排序
  head(zhnagxu2)
  colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)

  zhnagxu3 = zhnagxu2
  
  ##制作标签坐标，标签位于顶端
  # Taxonomies_x = ddply(zhnagxu3,"group", summarize, label_y = cumsum(Abundance))
  # head(Taxonomies_x )
  #标签位于中部
  # Taxonomies_x1 = ddply(zhnagxu3,"group", transform, label_y = cumsum(Abundance) - 0.5*Abundance)
  Taxonomies_x = plyr::ddply(zhnagxu3,"group", summarize,label_sd = cumsum(Abundance),label_y = cumsum(Abundance) - 0.5*Abundance)
  head( Taxonomies_x )
  
  # Taxonomies_x$label_y =
  Taxonomies_x = cbind(as.data.frame(zhnagxu3),as.data.frame(Taxonomies_x)[,-1])
  

  Taxonomies_x$label = Taxonomies_x$aa
  # #使用循环将堆叠柱状图柱子比较窄的别写标签，仅仅宽柱子写上标签
  # for(i in 1:nrow(Taxonomies_x)){
  #   if(Taxonomies_x[i,3] > 3){
  #     Taxonomies_x[i,5] = Taxonomies_x[i,5]
  #   }else{
  #     Taxonomies_x[i,5] = NA
  #   }
  # }

  
  
  Taxonomies_x$aa = factor(Taxonomies_x$aa,order = TRUE,levels = c(as.character(cc$aa)))
  
  
  
  ##普通柱状图
  
  
  p4 <- ggplot(Taxonomies_x , aes(x =  group, y = Abundance, fill = aa, order = aa)) +
    geom_bar(stat = "identity",width = 0.5,color = "black") +
    # scale_fill_manual(values = colors) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text=element_text(size=6)) +
    scale_y_continuous(name = "Relative abundance (%)") +
    guides(fill = guide_legend(title = j)) +
    labs(x="",y="Relative abundance (%)",
         title= "")
  # paste("Top ",Top," ",j,sep = "")
  p4
  if (is.na(axis_order)) {
    p4 = p4
  }else{
    p4 = p4 + scale_x_discrete(limits = axis_order)
  }
  
  
  if (sd == TRUE) {
    p4 =  p4 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
  }
  
  if (label == TRUE) {
    p4 = p4 +
      geom_text(aes(y = label_y, label = label ))
  }
  
  
  map = as.data.frame(phyloseq::sample_data(ps))
  if (length(unique(map$Group))>3){p4 = p4 + theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}
  
  #-------------冲击图
  cs = Taxonomies_x $aa
  
  lengthfactor <- cs %>%
    levels() %>%
    length()
  cs4 <- cs %>%
    as.factor() %>%
    summary()  %>%
    as.data.frame()
  cs4$id = row.names(cs4)
  
  
  #对因子进行排序
  df_arrange<- dplyr::arrange(cs4, id)
  #对Taxonomies_x 对应的列进行排序
  Taxonomies_x1<- dplyr::arrange(Taxonomies_x , aa)
  head(Taxonomies_x1)
  #构建flow的映射列Taxonomies_x
  Taxonomies_x1$ID = factor(rep(c(1:lengthfactor), cs4$.))
  
  #colour = "black",size = 2,,aes(color = "black",size = 0.8)
  head(Taxonomies_x1)
  Taxonomies_x1$Abundance
  
  p3 <- ggplot(Taxonomies_x1, aes(x = group, y = Abundance,fill = aa,alluvium = aa,stratum = ID)) +
    ggalluvial::geom_flow(aes(fill = aa, colour = aa),
              stat = "alluvium", lode.guidance = "rightleft",
              color = "black",size = 0.2,width = 0.35,alpha = .2)  +
    geom_bar(width = 0.45,stat = "identity") +
    labs(x="",y="Relative abundance (%)",
         title= "") +
    guides(fill = guide_legend(title = j),color = FALSE) +
    scale_y_continuous(expand = c(0,0))
  p3
  
  
  # flower plot
  p1 <- ggplot(Taxonomies_x1,
               aes(x = group, alluvium = aa, y = Abundance)) +
    ggalluvial::geom_flow(aes(fill = aa, colour = aa), width = 0)  +
    labs(x="",y="Relative abundance (%)",
         title="") +
    guides(fill = guide_legend(title = j),color = FALSE) +
    scale_y_continuous(expand = c(0,0))
  
  
  if (is.na(axis_order)) {
    p1 = p1
    p3 = p3
    
  }else{
    p1 = p1 + scale_x_discrete(limits = axis_order)
    p3 = p3 + scale_x_discrete(limits = axis_order)
    
  }
  # p3
  if (label == TRUE) {
    p1 = p1 +
      geom_text(aes(y = label_y, label = label ))
    p3 = p3 +
      geom_text(aes(y = label_y, label = label ))
  }
  
  if (sd == TRUE) {
    p1 =  p1 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
    p3 =  p3 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
  }
  
  if (length(unique(map$Group))>3){	p3=p3+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}
  
  return(list(p4,Taxonomies,p3,p1))
}



# 
# tax_glom_wt <- function(ps = ps,ranks = "Phylum") {
#   
#   
#   otu <- as.data.frame(t(vegan_otu(ps)))
#   tax <- as.data.frame(vegan_tax(ps))
#   
#   # building group
#   tax[[ranks]][is.na(tax[[ranks]])] = "Unknown"
#   tax[[ranks]][tax[[ranks]] == ""] = "Unknown"
#   split <- split(otu,tax[[ranks]])
#   #calculate sum by group
#   apply <- lapply(split,function(x)colSums(x[]))
#   # result chack
#   otucon <- do.call(rbind,apply)
#   
#   taxcon <- tax[1:match(ranks,colnames(tax))]
#   taxcon <- taxcon[!duplicated(tax[[ranks]]),]
#   #-tax name with NA wound be repeated with unknown
#   taxcon[[ranks]][is.na(taxcon[[ranks]])] = "unknown"
#   row.names(taxcon) <- taxcon[[ranks]]
#   
#   
#   pscon <- phyloseq(
#     otu_table( as.matrix(otucon),taxa_are_rows = TRUE),
#     tax_table(as.matrix(taxcon)),
#     sample_data(ps)
#   )
#   return(pscon)
# }
# 
