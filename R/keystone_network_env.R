


#
# res = keystone.micro(igraph, method = "zipi",select = "igraph.degree",  num= 10)

keystone.micro = function(
    igraph,
    method = "zipi",#"nodeproperty","huoscore","zipi"
    select = "igraph.degree",# "igraph.degree","igraph.closeness","igraph.betweenness","igraph.cen.degree"
    num= 10
){
  if (method == "zipi") {
    res = ZiPiPlot(igraph = igraph)
    p <- res[[1]]
    #p
    tem <- res[[2]] %>% dplyr::filter(roles != "Peripherals")
    id = row.names(tem)
    id
    dat1 = res[[2]]

  }else if(method == "nodeproperty"){
    #--使用节点属性判断

    node.p = node_properties(igraph)
    head(node.p)
    colnames(node.p)
    hub = node.p[,select] %>%
      sort(decreasing = TRUE) %>%
      head(num) %>%
      as.data.frame()
    colnames(hub) = "hub_sca"
    # head(hub)
    id = row.names(hub)

    dat1 = data.frame(row.names = id,id = id,roles = "Keystone")

  }else if(method == "huoscore"){

    #使用hubscore判断
    hub = hub_score(igraph)$vector %>%
      sort(decreasing = TRUE) %>%
      head(num) %>%
      as.data.frame()
    colnames(hub) = "hub_sca"
    head(hub)
    id = row.names(hub)
    dat1 = data.frame(row.names = id,id = id,roles = "Keystone")

  }

  return(list(id,dat1,p))
}






net_properties.4 <-function(igraph, n.hub = FALSE
){
  # igraph.weight <- E(igraph)$weight
  # network property
  # The size of the graph (number of edges)
  num.edges <- length(igraph::E(igraph))
  num.edges
  #  Order (number of vertices) of a graph
  num.vertices <- length(igraph::V(igraph))# length(diversity(igraph, weights = NULL, vids = 	V(igraph)))
  num.vertices

  connectance <- edge_density(igraph,loops=FALSE)

  # (Average degree)
  average.degree <- mean(igraph::degree(igraph))# 或者为2M/N,其中M 和N 分别表示网络的边数和节点数。
  average.degree
  # (Average path length)
  if (!is.null(igraph::E(igraph)$weight)) {
    igraph.weight <- igraph::E(igraph)$weight
    igraph::E(igraph)$weight = abs(igraph::E(igraph)$weight)
  }

  average.path.length <- average.path.length(igraph) # 同mean_distance(igraph) # mean_distance calculates the average path length in a graph
  average.path.length

  # (Diameter)
  diameter <- diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
  diameter


  if (!is.null(igraph::E(igraph)$weight)) {
    igraph::E(igraph)$weight = igraph.weight
  }
  #  edge connectivity / group adhesion
  edge.connectivity <- edge_connectivity(igraph)
  edge.connectivity
  # (Clustering coefficient)
  clustering.coefficient <- transitivity(igraph,type = "average")
  clustering.coefficient

  no.clusters <- no.clusters(igraph)
  no.clusters
  # (Degree centralization)
  centralization.degree <- centralization.degree(igraph)$centralization
  centralization.degree
  # (Betweenness centralization)
  centralization.betweenness <- centralization.betweenness(igraph)$centralization
  centralization.betweenness
  # (Closeness centralization)
  centralization.closeness <- centralization.closeness(igraph)$centralization
  centralization.closeness
  if (!is.null(igraph::E(igraph)$weight)) {
    num.pos.edges<-sum(igraph.weight>0)# number of postive correlation
    num.neg.edges<-sum(igraph.weight<0)# number of negative correlation
  }else{
    num.pos.edges<-0# number of postive correlation
    num.neg.edges<-0# number of negative correlation
  }

  #-----add RM
  modularity_igraph = function(net,method = "cluster_walktrap"){
    if (method == "cluster_walktrap" ) {
      fc <- cluster_walktrap(net)# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
    }

    if (method == "cluster_edge_betweenness" ) {
      fc <- cluster_edge_betweenness(net)# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
    }
    if (method == "cluster_fast_greedy" ) {
      fc <- cluster_fast_greedy(net)# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
    }
    if (method == "cluster_spinglass" ) {
      fc <- cluster_spinglass(net)# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
    }
    modularity <- modularity(net,membership(fc))
    return(modularity)
  }

  mod1 = modularity_igraph(igraph,method = "cluster_walktrap")
  rand.g <- erdos.renyi.game(length(igraph::V(igraph)), length(igraph::E(igraph)),type = "gnm")
  mod2 = modularity_igraph(rand.g,method = "cluster_walktrap")

  RM = (mod1-mod2)/mod2
  #---
  if (n.hub) {
    res = ZiPiPlot(igraph = igraph,method = "cluster_walktrap")
    data = res[[2]]
    head(data)
    n.hub = data$roles[data$roles != "Peripherals" ] %>% length()
  } else {
    n.hub = "Not.calculated"
  }


  igraph.network.pro <- rbind(num.edges,num.pos.edges,
                              num.neg.edges,num.vertices,
                              connectance,average.degree,average.path.length,
                              diameter,edge.connectivity,clustering.coefficient,
                              no.clusters,
                              centralization.degree,centralization.betweenness,
                              centralization.closeness,
                              RM,
                              mod1,
                              mod2,
                              n.hub

  )
  rownames(igraph.network.pro)<-c("num.edges(L)",
                                  "num.pos.edges",
                                  "num.neg.edges",
                                  "num.vertices(n)",
                                  "Connectance(edge_density)","average.degree(Average K)","average.path.length",
                                  "diameter","edge.connectivity",
                                  "mean.clustering.coefficient(Average.CC)",
                                  "no.clusters","centralization.degree",
                                  "centralization.betweenness","centralization.closeness",
                                  "RM(relative.modularity)",
                                  "modularity.net",
                                  "modularity_random",
                                  "the.number.of.keystone.nodes"
  )
  colnames(igraph.network.pro)<- "value"
  return(igraph.network.pro)
}




network.3 = function(
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
    keystone = T,
    methodkyestone = "zipi",
    clu_method = "cluster_fast_greedy",
    step = 100,
    yourmem = theme_void(),
    ncol = 3,
    nrow = 1,
    R = 10,
    ncpus = 1,
    random.net = F,
    envRDA = envRDA
){

  #--imput data
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  if (scale ) {ps_rela  = scale_micro(ps = ps,method = "rela")} else {ps_rela <- ps}
  mapping = as.data.frame(sample_data(ps_rela))
  y = matrix(1,nrow = 18,ncol = length(unique(mapping$Group)))
  layouts = as.character(unique(mapping$Group))
  mapping$ID = row.names(mapping)
  plots = list()
  plots1 = list()
  aa = 1
  # layout = layouts[1]
  for (layout in layouts) {

    mapi <- mapping[mapping$Group ==  layout,]
    psi = phyloseq(otu_table(ps_rela),
                   tax_table(ps_rela),
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
    #-----culculate edge
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

    # plotname = paste(path,"/network",layout,"_cover.pdf",sep = "")
    # ggsave(plotname, pnet1, width = 11, height = 9)
    plots[[aa]] = pnet
    plots1[[aa]] = pnet1

    # nodepro = node_properties(igraph)
    if (keystone ) {

      res = keystone.micro(igraph = igraph,
                           method = methodkyestone,
                           select = "igraph.degree",
                           num= 10)
      id = res[[1]]
      dat1 = res[[2]]
      p = res[[3]]
      #----culculate zi pi
      # res = ZiPiPlot(igraph = igraph,method = clu_method)
      # p <- res[[1]]
      ggsave(paste(path,"/",layout,"_ZiPi.pdf",sep = ""),p,width = 12, height = 10)
      # ZiPi <- res[[2]]
      write.csv(dat1 ,paste(path,"/",layout,"ZiPi.csv",sep = ""),row.names = FALSE)
      # tem <- dat1 %>% dplyr::filter(roles != "Peripherals")
      # id = row.names(tem)

      if (length(id) == 0) {
        write.csv(tem,paste(path,"/",layout,"ZiPi_no_keystone.csv",sep = ""),row.names = FALSE)
      } else {
        #--关键微生物堆叠柱状图
        result = keystone.barMainplot(ps,id,"Phylum")
        p1 = result[[1]]
        p2 = result[[2]]

        filename = paste(path,layout,"imformation_keystone",".pdf",sep = "")
        ggsave(filename,p1,width = 12,height = 8)

        filename = paste(path,layout,"imformation_keystone_flow",".pdf",sep = "")
        ggsave(filename,p2,width = 12,height = 8)

      }




      if (!is.null(envRDA)) {

        if (length(id) == 0) {
          write.csv(tem,paste(path,"/",layout,"No_keystone_cor_env.csv",sep = ""),row.names = FALSE)

        } else {
          #---关键微生物和环境因子相关
          result = zipi.keystone.cor.env (ps,
                                          id = id,
                                          envRDA = envRDA
          )
          p1 = result[[1]]
          filename = paste(path,layout,"ggheatmap_Keystone_env.pdf",sep = "")
          ggsave(filename,p1,width = 10,height = dim(envRDA)[2]/3)
        }

        #---基于zipi分类的四大类OTU矩阵和环境因子们特尔检验#----
        result = hub_mantel.env(ps,dat = dat1,envRDA)
        p1 = result[[1]]
        p2 = result[[2]]
        tab = result[[3]]

        if (is.null(p1)) {
          filename = paste(path,layout,"no_Keystone.csv",sep = "")
          write.csv(tab,filename)
        } else {
          filename = paste(path,layout,"Keystone.class_Mantel.csv",sep = "")
          write.csv(tab,filename)
          filename = paste(path,layout,"ggheatmap_Keystone.class.mantel.env.pdf",sep = "")
          ggsave(filename,p1,width = 4,height = 12)
          filename = paste(path,layout,"ggbubble_Keystone.class.mantel.env.pdf",sep = "")
          ggsave(filename,p2,width = 4,height = 12)
        }


      }



    }
    netpro_result<- net_properties.4(igraph)
    colnames(netpro_result)<-layout


    if (random.net) {
      result = random_Net_compate(igraph = igraph, type = "gnm", step = 100, netName = layout)
      p1 = result[[1]]
      sum_net = result[[4]]

      plotname = paste(path,"/Power_law_distribution_",layout,".pdf",sep = "")
      ggsave(plotname, p1, width = 8, height =6)

      write.csv(sum_net,paste(path,"/",layout,"_net_VS_erdos_properties.csv",sep = ""),row.names = TRUE)

    }
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


# head(dat)

hub_mantel.env = function(
    ps,
    dat,
    envRDA

){

  A = dat $roles %>% unique()
  A

  # if (length(A) != 1) {
  B = list()
  for (i in 1:length(A)) {
    t.1 <- dat %>% filter(roles == A[i]) %>% row.names()
    otu = phyloseq::otu_table(ps %>% scale_micro())
    tax = phyloseq::tax_table(ps%>% scale_micro())
    data = otu[t.1,] %>% #t() %>%
      as.data.frame()
    head(data)


    B[[i]] = data
  }
  names(B) = A
  tabOTU1 = B

  envRDA = envRDA[colnames(otu),]

  #--- mantel test
  rep = ggClusterNet::MetalTast (env.dat = envRDA, tabOTU = tabOTU1,
                                 distance = "bray",method = "metal")
  repR = rep[c(-seq(from=1,to=dim(rep)[2],by=2)[-1])]
  repP = rep[seq(from=1,to=dim(rep)[2],by=2)]
  head( repR)
  head( repP)
  mat = cbind(repR,repP)
  head(mat)


  pcm = reshape2::melt(repR, id = c("Envs"))
  head(pcm)
  pcm2 = reshape2::melt(repP, id = c("Envs"))
  head(pcm2)
  colnames(pcm2)[3] = "p"

  pcm2$lab = pcm2$p
  pcm2$lab[pcm2$lab < 0.001] = "**"
  pcm2$lab[pcm2$lab < 0.05] = "*"
  pcm2$lab[pcm2$lab >= 0.05] = ""
  pcm2$variable = NULL
  # pcm$variable = as.character(pcm$variable )
  # pcm2$variable = as.character(pcm2$variable)

  pcm3 = pcm %>% dplyr::left_join(pcm2)
  head(pcm3)
  colnames(pcm3)[1] = "id"
  # pcm3$p

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
      axis.text.x = element_text(colour = "black",angle = 90,vjust = 1,hjust = 1)
    )

  p1
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
      axis.text.x = element_text(colour = "black",angle = 90,vjust = 1,hjust = 1)
    )

  # } else{
  #
  #   p1 = NULL
  #   p2 = NULL
  #   pcm3 = NULL
  #
  # }

  return(list(p1,p2,pcm3))
}







# result = zipi.keystone.cor.env (ps,
#                                 id = id,
#                                 envRDA = env1
#                                 )
#
# result[[1]]

zipi.keystone.cor.env = function(
    ps,
    id,
    envRDA){

  otu = ps %>% vegan_otu() %>% t() %>%
    as.data.frame()
  ps.t = ps %>% scale_micro()
  otu = phyloseq::otu_table(ps.t)
  tax = phyloseq::tax_table(ps.t)
  head(otu)
  data = otu[id,] %>% t() %>%
    as.data.frame()
  match(row.names(data),row.names(envRDA))
  envRDA = envRDA[sample_names(ps),]

  result = cor_env_ggcorplot(
    env1 = envRDA,
    env2 = data,
    label =  F,
    col_cluster = F,
    row_cluster = F,
    method = "spearman",
    r.threshold= 0,
    p.threshold= 0

  )

  p1 <- result[[1]]
  p1
  p2 <- result[[2]]
  p2

  p0 = p1|p2
  # filename = paste(path,layout,"Keystone_abundance.csv",sep = "")
  # write.csv(data,filename)
  #
  # filename = paste(path,layout,"ggheatmap_Keystone_env.pdf",sep = "")
  # ggsave(filename,p1,width = 4,height = dim(envRDA)[2]/3)
  # filename = paste(path,layout,"ggbubble_Keystone_env.pdf",sep = "")
  # ggsave(filename,p2,width = 4,height = dim(envRDA)[2]/3)

  return(list(p0))
}















# tem <- res[[2]] %>% dplyr::filter(roles != "Peripherals")
# id = row.names(tem)
# id
# result = keystone.barMainplot(ps,id,"Phylum")
# result[[1]]
# result[[2]]


keystone.barMainplot = function(
    ps,
    id,
    j = "Family"){
  otu = ps %>% scale_micro() %>%
    vegan_otu() %>% t() %>%
    as.data.frame()
  ps.t = ps %>% scale_micro()
  otu_table(ps.t) = otu_table(as.matrix(otu[id,]),taxa_are_rows = TRUE)
  #--关键微生物和土壤养分关系

  #--丰度展示
  # source("E:\\Shared_Folder\\Function_local\\R_function\\micro/barMainplot.R")

  result = barMainplot(ps = ps.t,
                       tran = F,
                       j = j,
                       # axis_ord = axis_order,
                       label = FALSE,
                       sd = FALSE,
                       Top = 10)
  p4_1 <- result[[1]] +
    scale_fill_hue() + theme_classic()
  p4_1

  # filename =  paste(path,layout,"imformation_keystone",".pdf",sep = "")
  # ggsave(filename,p4_1,width = 6,height = 8)

  p4_2  <- result[[3]]  +
    scale_fill_hue() + theme_classic()
  p4_2
  # filename = paste(path,layout,"imformation_keystone_flow",".pdf",sep = "")
  # ggsave(filename,p4_2,width = 6,height = 8)

  result = barMainplot(ps = ps.t,
                       tran = T,
                       j = j,
                       # axis_ord = axis_order,
                       label = FALSE,
                       sd = FALSE,
                       Top = 10)
  p3_1 <- result[[1]] +
    scale_fill_hue() + theme_classic()
  # p3_1
  # filename = paste(path,layout,"imformation_keystone_100",".pdf",sep = "")
  # ggsave(filename,p3_1,width = 6,height = 8)

  p3_2  <- result[[3]]  +
    scale_fill_hue() + theme_classic()
  # p3_2
  p1 = p4_1|p3_1
  p2 = p4_2|p3_2
  # filename = paste(path,layout,"imformation_keystone_100_folw",".pdf",sep = "")
  # ggsave(filename,p3_2,width = 6,height = 8)

  return(list(p1,p2,ps.t))
}
