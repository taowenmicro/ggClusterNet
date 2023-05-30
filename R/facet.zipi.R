
#facet.zipi
facet.zipi = function(
    net.dat.4= net.dat.4
){
  x1<- c(0, 0.62,0,0.62)
  x2<- c( 0.62,1,0.62,1)
  y1<- c(-Inf,2.5,2.5,-Inf)
  y2 <- c(2.5,Inf,Inf,2.5)
  lab <- c("peripheral",'Network hubs','Module hubs','Connectors')

  roles.colors <- c("#E6E6FA","#DCDCDC","#F5FFFA", "#FAEBD7")


  tab = data.frame(x1 = x1,y1 = y1,x2 = x2,y2 = y2,lab = lab)
  tem = net.dat.4$group %>% unique() %>% length()
  for ( i in 1:tem) {
    if (i == 1) {
      tab2 = tab
    } else{
      tab2 = rbind(tab2,tab)
    }
  }


  p <- ggplot() +
    geom_rect(data=tab2,
              mapping=aes(xmin=x1,
                          xmax=x2,
                          ymin=y1,
                          ymax=y2,
                          fill = lab))+
    guides(fill=guide_legend(title="Topological roles")) +
    scale_fill_manual(values = roles.colors)+
    geom_point(data=net.dat.4,aes(x=p, y=z,color=module)) + theme_bw()+
    guides(color= F) +
    ggrepel::geom_text_repel(data = net.dat.4,
                             aes(x = p, y = z,
                                 color = module,label=label),size=4)+
    # facet_wrap(.~group) +
    facet_grid(.~ group, scale='free') +
    theme(strip.background = element_rect(fill = "white"))+
    xlab("Participation Coefficient")+ylab(" Within-module connectivity z-score")
  return(p)
}
