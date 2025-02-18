
#' Microbial related network
#'
#' @param ps phyloseq Object, contains OTU tables, tax table and map table, represented sequences,phylogenetic tree.
#' @param N filter OTU tables by abundance.The defult, N=0, extract the top N number relative abundance of OTU. e.g 100
#' @param r.threshold The defult, r.threshold=0.6, it represents the correlation that the absolute value
#'  of the correlation threshold is greater than 0.6. the value range of correlation threshold from 0 to 1.
#' @param p.threshold The defult, p.threshold=0.05, it represents significance threshold below 0.05.
#' @param select_layout  TURE or FALSE
#' @param layout_net defult "model_maptree"
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
#' @param layout_net select layout from ggClusterNet
#' @param big TRUE or FALSE the number of micro data was so many (> 300),you can chose TREU
#' @examples
#' data(ps)
#' path = "./netowrk/"
#' dir.create(path)
#' result = network.2(ps = ps,N = 100,r.threshold=0.6,big = T,
#'                    select_layout = T,
#'                    p.threshold=0.05,label = FALSE,path = path ,zipi = F)
#' result[[1]]
#' result[[2]]
#' @return list which contains OTU correlation matrix
#' @author Contact: Tao Wen \email{taowen@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn} yongxin liu \email{yxliu@@genetics.ac.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Tao Wen#, Penghao Xie#, Shengdie Yang, Guoqing Niu, Xiaoyu Liu, Zhexu Ding, Chao Xue, Yong-Xin Liu *, Qirong Shen, Jun Yuan*
#' ggClusterNet: an R package for microbiome network analysis and modularity-based multiple network layouts
#' iMeta 2022,DOI: \url{doi: 10.1002/imt2.32}
#' @export




network.2 = function(
  otu = NULL,
                   tax = NULL,
                   map = NULL,
                   ps = NULL,
                   N = 0,
                   big = FALSE,
                   select_layout = FALSE,
                   layout_net = "model_maptree2",
                   r.threshold = 0.6,
                   p.threshold = 0.05,
                    maxnode = 2,
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
  if (scale ) {ps_rela  = scale_micro(ps = ps,method = "rela")} else {ps_rela <- ps}
  mapping = as.data.frame(sample_data(ps_rela))
  y = matrix(1,nrow = 16,ncol = length(unique(mapping$Group)))
  layouts = as.character(unique(mapping$Group))
  mapping$ID = row.names(mapping)
  plots = list()
  plots1 = list()
  aa = 1
  # layout = layouts[1]
  for (layout in layouts) {

    mapi <- mapping[mapping$Group ==  layout,]
    psi = phyloseq(otu_table(ps_rela),
                   phyloseq::tax_table(ps_rela),
                   sample_data(mapi)
    ) %>%
      filter_OTU_ps(Top = N) %>%
      filter_taxa( function(x) sum(x ) > 0 , TRUE)
    print(layout)
    if (big == TRUE) {
      result = cor_Big_micro(ps = psi,N = 0,r.threshold= r.threshold,p.threshold=p.threshold,method = method,scale = FALSE)
      a = 2} else if(big == FALSE){
      result = corMicro (ps = psi,N = 0,r.threshold= r.threshold,p.threshold=p.threshold,method = method,R = R,ncpus = ncpus)
      a = 1}

    print("cor matrix culculating over")
    cor = result[[1]]    #Extract correlation matrix


    if (cor %>% as.vector() %>% max() == 0) {
      stop("The connect value in cor matrix all was zone")
    }


    if (select_layout) {
      node = NULL
      node = culculate_node_axis(
        cor.matrix = cor,
        layout = layout_net,
        seed = 1,
        group = NULL,
        model = FALSE,
        method = clu_method)

    }else {
      result2 <- model_Gephi.2(cor = cor,
                               method = clu_method,
                               seed = 12
      )
      node = result2[[1]]
    }

    # ps_net = psi
    otu_table = as.data.frame(t(vegan_otu(psi)))
    tax_table = as.data.frame(vegan_tax(psi))
    nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
    #-----culculate edge #--------
    edge = edgeBuild(cor = cor,node = node)
    #-------output---edges and nodes--to Gephi --imput--
    edge_Gephi = data.frame(source = edge$OTU_1,target = edge$OTU_2,correlation =  edge$weight,direct= "undirected",cor =  edge$cor)
    # building node table
    node_Gephi = data.frame(ID= nodes$elements,nodes[4:dim(nodes)[2]],Label = nodes$elements)

    idedge <- c(as.character(edge_Gephi$source),as.character(edge_Gephi$target))
    idedge <- unique(idedge)
    row.names(node_Gephi) <- as.character(node_Gephi$ID)
    node_Gephi1 <- node_Gephi[idedge, ]

    write.csv(edge_Gephi ,paste(path,"/",layout,"_Gephi_edge.csv",sep = ""),row.names = FALSE,quote = FALSE)
    write.csv(node_Gephi,paste(path,"/",layout,"_Gephi_allnode.csv",sep = ""),row.names = FALSE,quote = FALSE)
    write.csv(node_Gephi1,paste(path,"/",layout,"_Gephi_edgenode.csv",sep = ""),row.names = FALSE,quote = FALSE)


    igraph  = igraph::graph_from_data_frame(nodeEdge(cor = cor)[[1]], directed = FALSE, vertices = nodeEdge(cor = cor)[[2]])
    nodepro = node_properties(igraph)
    write.csv(nodepro,paste(path,"/",layout,"_node_properties.csv",sep = ""),row.names = TRUE)

    nodeG = merge(nodes,nodepro,by = "row.names",all.x  = TRUE)
    row.names(nodeG) = nodeG$Row.names
    nodeG$Row.names = NULL

    numna = (dim(nodeG)[2] - 3) : dim(nodeG)[2]
    nodeG[,numna][is.na(nodeG[,numna])] = 0

    # print(dim(nodeG))
    head(nodeG)
    # dat.alpha = nodeG[nodeG$igraph.degree == 0,]

    pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor),
                                    data = edge, size = 0.03,alpha = 0.5) +
      geom_point(aes(X1, X2,fill = !!sym(fill),size = !!sym(size) ),
                 pch = 21, data = nodeG,color = "gray40") +
      labs( title = paste(layout,"network",sep = "_"))  +
      # geom_text_repel(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
      # geom_text(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
      scale_colour_manual(values = c("#6D98B5","#D48852")) +
      scale_size(range = c(0.8, maxnode)) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
      theme(panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5)
            ) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      theme(legend.background = element_rect(colour = NA)) +
      theme(panel.background = element_rect(fill = "white",  colour = NA)) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())


    pnet1 <- ggplot() + geom_curve(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor),
                                   data = edge, size = 0.03,alpha = 0.3,curvature = -0.2) +
      geom_point(aes(X1, X2,fill = !!sym(fill),size = !!sym(size)),
                 pch = 21, data = nodeG,color = "gray40") +
      labs( title = paste(layout,"network",sep = "_"))  +
      # geom_text_repel(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
      # geom_text(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
      scale_colour_manual(values = c("#6D98B5","#D48852")) +
      scale_size(range = c(0.8,maxnode)) +
      scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = NULL) +
      theme(panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      theme(legend.background = element_rect(colour = NA)) +
      theme(panel.background = element_rect(fill = "white",  colour = NA)) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
    pnet1



    if (label == TRUE ) {
      pnet <- pnet +  ggrepel::geom_text_repel(aes(X1, X2,label=!!sym(lab)),size=4, data = nodeG)
      pnet1 <- pnet1 +  ggrepel::geom_text_repel(aes(X1, X2,label=!!sym(lab)),size=4, data = nodeG)
    }

    plotname = paste(path,"/network",layout,".pdf",sep = "")
    ggsave(plotname, pnet, width = 11, height = 9)

    plotname = paste(path,"/network",layout,"_cover.pdf",sep = "")
    ggsave(plotname, pnet1, width = 11, height = 9)
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
    netpro_result<- net_properties.2(igraph)
    colnames(netpro_result)<-layout


    result = random_Net_compate(igraph = igraph, type = "gnm", step = 100, netName = layout)
    p1 = result[[1]]
    sum_net = result[[4]]

    plotname = paste(path,"/Power_law_distribution_",layout,".pdf",sep = "")
    ggsave(plotname, p1, width = 8, height =6)

    write.csv(sum_net,paste(path,"/",layout,"_net_VS_erdos_properties.csv",sep = ""),row.names = TRUE)

    y = as.data.frame(y)
    colnames(y) = layouts
    # head(y)
    y[layout] = netpro_result[,1]
    row.names(y) = row.names(netpro_result)
    aa = aa+1
  }

  plotname = paste(path,"/network_all.pdf",sep = "")
  p  = ggpubr::ggarrange(plotlist = plots,
                         common.legend = TRUE, legend="right",ncol = ncol,nrow = nrow)
  p1  = ggpubr::ggarrange(plotlist = plots1,
                          common.legend = TRUE, legend="right",ncol = ncol,nrow = nrow)

  if (length(layouts) == 1) {
    p = pnet
    p1 = pnet1
  }
  return(list(p,y,p1,cor))
}






# data(igraph)
# result = random_Net_compate(igraph = igraph, type = "gnm", step = 100, netName = "KO")
# p1 = result[[1]]
# sum_net = result[[4]]


random_Net_compate = function(
  igraph = igraph,
  type = "gnm",
  step = 100,
  netName = "KO"
){
  #--
  res = random_Net(igraph = igraph,type = "gnm")
  p = res[[1]]

  data = grobal_pro_compare(igraph = igraph,type = "gnm", step = 100, netName = "KO")

  return(list(p,plotdata = res[[2]],randomigraph = res[[3]],pro_compare = data ))
}


#--make random network compared with actral network
random_Net = function(
  igraph = igraph,
  type = "gnm"
){
  rand.g<- erdos.renyi.game(length(V(igraph)), length(E(igraph)),type = c(type))
  # degree_distribution

  ## We compare the degree of random network and this network
  data1 = data.frame(network= degree_distribution(igraph, cumulative = FALSE),group = "network",ID = c(1:length(degree_distribution(igraph, cumulative = FALSE))))
  data2 = data.frame(network = degree_distribution(rand.g, cumulative = FALSE) ,group = "E–R network",ID = c(1:length(degree_distribution(rand.g, cumulative = FALSE) )))
  data = rbind(data1,data2)
  p1 <- ggplot(data) +geom_point(aes(x = ID,y = network,group =group,fill = group),pch = 21,size = 2) +
    geom_smooth(aes(x = ID,y = network,group =group,color = group))+
    theme_bw() + theme(
      plot.margin=unit(c(0,0,0,0), "cm")
    )
  return(list(p1,plotdata = data,randomigraph = rand.g))
}

# plotname = paste(path,"/Power_law_distribution_",layout,".pdf",sep = "")
# ggsave(plotname, p1, width = 8, height =6)



grobal_pro_compare = function(
  igraph = igraph,
  type = "gnm",
  step = 100,
  netName = "KO"

){

  rand.g.netpro_result<-c()

  for (i in 1:step){
    #####random null model
    rand.g <- erdos.renyi.game(length(igraph::V(igraph)), length(igraph::E(igraph)),type = c(type))
    tem_netpro_result<-net_properties.2.rm(rand.g)
    tem_netpro_result[16,1] = 0
    tem_netpro_result = as.data.frame(tem_netpro_result)
    tem_netpro_result$value = as.numeric(tem_netpro_result$value)
    tem_netpro_result = as.matrix(tem_netpro_result)
    rand.g.netpro_result<-cbind(rand.g.netpro_result,tem_netpro_result)
  }


  result_summary<-cbind(rowMeans(rand.g.netpro_result),apply(rand.g.netpro_result,1,sd))
  colnames(result_summary)<-c("Means","SD")
  head(result_summary)
  netpro_result<- net_properties.2(igraph)
  colnames(netpro_result)<- netName
  sum_net = cbind(netpro_result,result_summary)
  return(sum_net)
}











