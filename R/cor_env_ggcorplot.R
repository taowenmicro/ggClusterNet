



cor_env_ggcorplot <- function(
  env1 = env1,
  env2 = env2,
  label =  T,
  col_cluster = T,
  row_cluster = T,
  method = "spearman",
  r.threshold=0.6,
  p.threshold=0.05,
  theme.size = 10
  
){

  if (dim(env2)[2] == 1) {
    env2 = env2
  } else {
    env2 <- env2[match(row.names(env1),row.names(env2)),]
  }

  env0 <- cbind(env1,env2)
  occor = psych::corr.test(env0,use="pairwise",method=method,adjust="fdr",alpha=.05)
  occor.r = occor$r
  occor.p = occor$p
  
  occor.r[occor.p > p.threshold&abs(occor.r) < r.threshold] = 0
  
  head(env0)
  # data[data > 0.3]<-0.3
  #drop gene column as now in rows
  
  if (col_cluster) {
    clust <- hclust(dist(env1 %>% as.matrix()%>% t())) # hclust with distance matrix
    ggtree_plot <- ggtree::ggtree(clust)
  }
  if (row_cluster) {
    v_clust <- hclust(dist(env2 %>% as.matrix() %>% t()))
    ggtree_plot_col <- ggtree::ggtree(v_clust) + ggtree::layout_dendrogram()
  }
  
  
  occor.r = as.data.frame(occor.r)
  
  if (dim(env2)[2] == 1) {
    
    data <- occor.r[colnames(env1),colnames(env2)]
    data = data.frame(row.names = colnames(env1),data)
    colnames(data) = colnames(env2)
    data$id = row.names(data)
  } else {
    data <- occor.r[colnames(env1),colnames(env2)]
    data$id = row.names(data)
  }
  
  

  pcm = reshape2::melt(data, id = c("id"))
  head(pcm)
  occor.p = as.data.frame(occor.p)
  
  if (dim(env2)[2] == 1) {
    
    data <- occor.p[colnames(env1),colnames(env2)]
    data = data.frame(row.names = colnames(env1),data)
    colnames(data) = colnames(env2)
    data$id = row.names(data)
  } else {
    data <- occor.p[colnames(env1),colnames(env2)]
    data$id = row.names(data)
    
  }
  

  pcm2 = reshape2::melt(data, id = c("id"))
  head(pcm2)
  colnames(pcm2)[3] = "p"
  
  pcm2$lab = pcm2$p 
  pcm2$lab[pcm2$lab < 0.001] = "**"
  pcm2$lab[pcm2$lab < 0.05] = "*"
  pcm2$lab[pcm2$lab >= 0.05] = ""
  pcm3 = pcm %>% left_join(pcm2)
  
  
  p1 = ggplot(pcm3, aes(y = id, x = variable)) + 
    # geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) + 
    geom_tile(aes(size = value,fill = value))+
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) + 
    geom_text(aes(label = lab)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  + 
    # scale_fill_manual(values = colours, guide = FALSE) + 
    scale_x_discrete(limits = rev(levels(pcm$variable)))  + 
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colours =colorRampPalette(c("#377EB8","#F7F4F9","#E41A1C"))(60)) +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",size = theme.size,angle = 60,vjust = 1,hjust = 1)
    )
  

  p2 = ggplot(pcm3, aes(y = id, x = variable)) + 
    geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) + 
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) + 
    geom_text(aes(label = lab)) +
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  + 
    # scale_fill_manual(values = colours, guide = FALSE) + 
    scale_x_discrete(limits = rev(levels(pcm$variable)))  + 
    scale_y_discrete(position = "right")  +
    scale_fill_gradientn(colours =colorRampPalette(c("#377EB8","#F7F4F9","#E41A1C"))(60))  +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",size = theme.size,angle = 60,vjust = 1,hjust = 1)
    )
  
  if (col_cluster) {
    p1 <- p1  %>%
      aplot::insert_left(ggtree_plot, width=.2) 
    p2 <- p2  %>%
      aplot::insert_left(ggtree_plot, width=.2) 
  }
  
  if (label) {
    p1 <- p1  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
    p2 <- p2  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
  }
  return(list(p1,p2))
}
