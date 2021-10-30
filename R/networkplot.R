
#' Microbial related network
#'
#' @param ps phyloseq Object, contains OTU tables, tax table and map table, represented sequences,phylogenetic tree.
#' @param N filter OTU tables by abundance.The defult, N=0, extract the top N number relative abundance of OTU.
#' @param r.threshold The defult, r.threshold=0.6, it represents the correlation that the absolute value
#'  of the correlation threshold is greater than 0.6. the value range of correlation threshold from 0 to 1.
#' @param p.threshold The defult, p.threshold=0.05, it represents significance threshold below 0.05.
#' @param method method for Correlation calculation,method="pearson" is the default value. The alternatives to be passed to cor are "spearman" and "kendall".
#' @param label Whether to add node label.
#' @param group Separate Group.
#' @param lay layout which network show
#' @param path save path of all of network analyse.
#' @param fill fill coulor of node
#' @param size node size
#' @param zipi zipi Calculation
#' @param step Random network sampling times
#' @param R repeat number of p value calculate
#' @param ncpus number of cpus used for sparcc
#' @examples
#' data(ps)
#' result = network (ps = ps,N = 0.001,r.threshold=0.6,p.threshold=0.05,label = FALSE,path = path ,zipi = TRUE)
#' result[[1]]
#' result[[2]]
#' @return list which contains OTU correlation matrix
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export




network = function(otu = NULL,
                   tax = NULL,
                   map = NULL,
                   ps = NULL,
                   N = 0,
                   r.threshold = 0.6,
                   p.threshold = 0.05,
                   method = "spearman",
                   label = FALSE,
                   lab = "elements",
                   group = "Group",
                   path = "./",
                   fill = "Phylum",
                   size = "igraph.degree",
                   scale = TRUE,
                   zipi = FALSE,
                   clu_method = "cluster_fast_greedy",
                   step = 100,
                   yourmem = theme_void(),
                   ncol = 3,
                   nrow = 1,
                   R = 10,
                   ncpus = 1
                   ){


  #--imput data ---------
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  #---------------------------data washing-------------------------------------------------
  #transform to relative abundance
  if (scale == TRUE) {
    ps_rela  = scale_micro(ps = ps,method = "rela")

  } else {
    ps_rela <- ps
  }

  mapping = as.data.frame(sample_data(ps))
  ps_sub = filter_OTU_ps(ps = ps,Top = N)
  y = matrix(1,nrow = 14,ncol = length(unique(mapping$Group)))
  #--transmit N
  # d = N

  layouts = as.character(unique(mapping$Group))
  mapping$ID = row.names(mapping)
  ##################---------------------------------------calculate network---------------------------------------------------

  plots = list()
  plots1 = list()
  # layout = layouts[1]
  aa = 1
  for (layout in layouts) {
    mapi <- mapping[mapping$Group ==  layout,]
    psi = phyloseq(otu_table(ps_sub),
                   tax_table(ps_sub),
                   sample_data(mapi )
                   )

    psi = filter_taxa(psi, function(x) sum(x ) > 0 , TRUE)


    print(layout)
    result = corMicro (ps = psi,N = 0,r.threshold= r.threshold,p.threshold=p.threshold,method = method,R = R,ncpus = ncpus)
    print("cor matrix culculating over")
    cor = result[[1]]    #Extract correlation matrix

    if (cor %>% as.vector() %>% max() == 0) {
      stop("The connect value in cor matrix all was zone")
    }


    result2 <- model_Gephi.2(cor = cor,
                             method = clu_method,
                             seed = 12
    )
    node = result2[[1]]
    dim(node)

    ps_net = result[[3]]
    otu_table = as.data.frame(t(vegan_otu(ps_net)))


    tax_table = as.data.frame(vegan_tax(ps_net))
    #---node节点注释#-----------
    nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
    head(nodes)

    #-----计算边#--------
    edge = edgeBuild(cor = cor,plotcord = node)
    head(edge)
    edge$weight
    colnames(edge)[8] = "cor"
    #-------output---edges and nodes--to Gephi --imput--
    edge_Gephi = data.frame(source = edge$OTU_1,target = edge$OTU_2,correlation =  edge$weight,direct= "undirected",cor =  edge$cor)
    # building node table
    node_Gephi = data.frame(ID= nodes$elements,nodes[4:dim(nodes)[2]],Label = nodes$elements)

    idedge <- c(as.character(edge_Gephi$source),as.character(edge_Gephi$target))
    idedge <- unique(idedge)
    row.names(node_Gephi) <- as.character(node_Gephi$ID)
    node_Gephi1 <- node_Gephi[idedge, ]

    write.csv(edge_Gephi ,paste(path,"/",layout,"_Gephi_edge.csv",sep = ""),row.names = FALSE)
    write.csv(node_Gephi,paste(path,"/",layout,"_Gephi_allnode.csv",sep = ""),row.names = FALSE)
    write.csv(node_Gephi1,paste(path,"/",layout,"_Gephi_edgenode.csv",sep = ""),row.names = FALSE)

    # a = nodeEdge(cor = cor)[[1]]
    # dim(a)
    # head(a)
    # a %>% filter(weight != 0)
    # as.vector(lower.tri(cor))
    igraph  = igraph::graph_from_data_frame(nodeEdge(cor = cor)[[1]], directed = FALSE, vertices = nodeEdge(cor = cor)[[2]])
    nodepro = node_properties(igraph)
    write.csv(nodepro,paste(path,"/",layout,"_node_properties.csv",sep = ""),row.names = TRUE)
    nodeG = merge(nodes,nodepro,by = "row.names",all.x  = TRUE)
    row.names(nodeG) = nodeG$Row.names
    nodeG$Row.names = NULL

    numna = (dim(nodeG)[2] - 3) : dim(nodeG)[2]
    nodeG[,numna][is.na(nodeG[,numna])] = 0

    head(nodeG)
    pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor),
                                    data = edge, size = 0.5,alpha = 0.3) +
      geom_point(aes(X1, X2,fill = !!sym(fill),size = !!sym(size)),pch = 21, data = nodeG) +
      labs( title = paste(layout,"network",sep = "_"))  +
      # geom_text_repel(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
      # geom_text(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
      scale_colour_manual(values = c("#377EB8","#E41A1C")) +
      scale_size(range = c(4, 14)) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
      theme(panel.background = element_blank()) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      theme(legend.background = element_rect(colour = NA)) +
      theme(panel.background = element_rect(fill = "white",  colour = NA)) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
    pnet

    pnet1 <- ggplot() + geom_curve(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor),
                                              data = edge, size = 0.5,alpha = 0.3,curvature = -0.2) +
      geom_point(aes(X1, X2,fill = !!sym(fill),size = !!sym(size)),pch = 21, data = nodeG) +
      labs( title = paste(layout,"network",sep = "_"))  +
      # geom_text_repel(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
      # geom_text(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
      scale_colour_manual(values = c("#377EB8","#E41A1C")) +
      scale_size(range = c(4, 14)) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
      theme(panel.background = element_blank()) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      theme(legend.background = element_rect(colour = NA)) +
      theme(panel.background = element_rect(fill = "white",  colour = NA)) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
    pnet1

    if (label == TRUE ) {
      pnet <- pnet +  geom_text_repel(aes(X1, X2,label=!!sym(lab)),size=4, data = nodeG)
      pnet1 <- pnet1 +  geom_text_repel(aes(X1, X2,label=!!sym(lab)),size=4, data = nodeG)
    }

    plotname = paste(path,"/network",layout,".pdf",sep = "")
    # ggsave(plotname, pnet, width = 16 * dim(nodeG)[1]/200, height = 14*dim(nodeG)[1]/200)
    ggsave(plotname, pnet, width = 16 , height = 14)

    plotname = paste(path,"/network",layout,"_cover.pdf",sep = "")
    # ggsave(plotname, pnet, width = 16 * dim(nodeG)[1]/200, height = 14*dim(nodeG)[1]/200)
    ggsave(plotname, pnet1, width = 16 , height = 14)
    plots[[aa]] = pnet
    plots1[[aa]] = pnet1

    # nodepro = node_properties(igraph)
    if (zipi ) {
      #----culculate zi pi
      res = ZiPiPlot(igraph = igraph,method = clu_method)
      p <- res[[1]]
      ggsave(paste(path,"/",layout,"_ZiPi.pdf",sep = ""),p,width = 12, height = 10)
      ZiPi <- res[[2]]
      write.csv(ZiPi ,paste(path,"/",layout,"ZiPi.csv",sep = ""),row.names = FALSE)
    }

    #----
    rand.g<- erdos.renyi.game(length(V(igraph)), length(E(igraph)),type = c("gnm"))
    # degree_distribution

    ## We compare the degree of random network and this network
    data1 = data.frame(network= degree_distribution(igraph, cumulative = FALSE),group = "E–R network",ID = c(1:length(degree_distribution(igraph, cumulative = FALSE))))
    data2 = data.frame(network = degree_distribution(rand.g, cumulative = FALSE) ,group = "network",ID = c(1:length(degree_distribution(rand.g, cumulative = FALSE) )))
    data = rbind(data1,data2)
    p1 <- ggplot(data) +geom_point(aes(x = ID,y = network,group =group,fill = group),pch = 21,size = 2) +
      geom_smooth(aes(x = ID,y = network,group =group,color = group))+
      theme_bw() + theme(
        plot.margin=unit(c(0,0,0,0), "cm")
      )
    plotname = paste(path,"/Power_law_distribution_",layout,".pdf",sep = "")
    ggsave(plotname, p1, width = 8, height =6)

    #---------------------

    rand.g.netpro_result<-c()
    for (i in 1:step){
      #####random null model
      rand.g<- erdos.renyi.game(length(V(igraph)), length(E(igraph)),type = c("gnm"))
      tem_netpro_result<-net_properties(rand.g)
      rand.g.netpro_result<-cbind(rand.g.netpro_result,tem_netpro_result)
    }
    #--------计算step次的结果
    # rand.g.netpro_result

    # --简化结果--------对随机矩阵结果求均值和sd值
    result_summary<-cbind(rowMeans(rand.g.netpro_result),apply(rand.g.netpro_result,1,sd))
    colnames(result_summary)<-c("Means","SD")
    head(result_summary)


    ###网络边的赋值及其设置
    igraph.weight <- E(igraph)$weight# 将igraph weight属性赋值到igraph.weight,用于后边做图
    # E(igraph)$weight <- NA
    # igraph<-remove.edge.attribute(igraph,"weight")#把边值删除
    netpro_result<- net_properties(igraph)
    colnames(netpro_result)<-layout

    sum_net = cbind(netpro_result,result_summary)
    write.csv(sum_net,paste(path,"/",layout,"_net_VS_erdos_properties.csv",sep = ""),row.names = TRUE)

    y = as.data.frame(y)
    colnames(y) = layouts
    # head(y)
    y[layout] = netpro_result[,1]
    row.names(y) = row.names(netpro_result)
    aa = aa+1
  }
  plotname = paste(path,"/network_all.pdf",sep = "")
  p  = ggpubr::ggarrange(plotlist = plots, common.legend = TRUE, legend="right",ncol = ncol,nrow = nrow)
  p1  = ggpubr::ggarrange(plotlist = plots1, common.legend = TRUE, legend="right",ncol = ncol,nrow = nrow)

  if (length(layouts) == 1) {
    p = pnet
    p1 = pnet1
  }
  return(list(p,y,p1))
}



