


ps.group.mean2 = function(ps,group = "Group",scale = TRUE){
  otu_table = as.data.frame(t(vegan_otu(ps)))
  head(otu_table)
  
  design = as.data.frame(sample_data(ps))
  ## 计算相对丰度，计算每个物种丰度均值，按照均值排序
  OTU = as.matrix(otu_table)
  if (scale == TRUE) {
    norm = t(t(OTU)/colSums(OTU,na=TRUE)) #* 100 # normalization to total 100
  } else{
    norm = OTU
  }
  norma = norm %>%
    t() %>% as.data.frame()
  #数据分组计算平均值
  iris.split <- split(norma,as.factor(design$Group))
  
  iris.apply <- lapply(iris.split,function(x)colMeans(x,na.rm = TRUE))
  # 组合结果
  norm2 <- do.call(rbind,iris.apply)%>%
    t()
  norm2 = as.data.frame(norm2)
  norm2$mean=apply(norm2,1,mean)
  norm2$ID = row.names(norm2)
  return(norm2)
}
