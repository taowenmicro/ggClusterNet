

#  ggClusterNet 2.0: an R package for microbial co-occurrence networks and associated indicator correlation patterns



![](https://github.com/taowenmicro/Rcoding/blob/main/Fig0.jpg?raw=true)


##  Introduce
Since the last version release in 2022, ggClusterNet has emerged as a critical resource for microbiome research, enabling microbial co-occurrence network analysis and visualization in over 200 studies (Google Scholar citations). To address emerging challenges in microbiome studies, including multi-factor experimental designs, multi-treatment, and multi-omics data, we present a comprehensive upgrade with the following four components: 1) We recommended and designed a microbial co-occurrence network analysis pipeline incorporating network computation and visualization (Pearson/Spearman/SparCC correlations), topological characterization of network and node properties, multi-network structure comparison and statistical testing, exploration of network stability (robustness), and identification and analysis of network modules; 2) Developed microbial network mining functions for multi-factor, multi-treatment, and spatiotemporal-scale analysis, such as Facet. Network(), module.compare.m.ts(), Robustness.Random.removal.ts(), etc.; 3) Developed functions for microbial and multi-factor interaction analysis, along with versatile visualization layout algorithms, such as MatCorPlot2(), Miccorplot3(), cor_link3(), matcorplotj(), and two.cor(); 4) Developed functions for cross-domain and multi-omics integrated network analysis, including corBionetwork.st(), and developed a comprehensive suite of visualization layout algorithms specifically designed for exploring complex relationships in these networks, such as model_maptree2(), model_Gephi.3(), cir.squ(), and cir.maptree2(). Collectively, the latest updates to ggClusterNet 2.0 empower researchers to explore complex network interactions with enhanced capabilities, offering a robust, efficient, user-friendly, reproducible, and visually versatile tool for microbial co-occurrence networks and associated indicator correlation patterns. The ggClusterNet 2.0 R package is open-source and freely accessible on GitHub (https://github.com/taowenmicro/ggClusterNet).


## Install

```


install.packages("BiocManager")
library(BiocManager)
install("remotes")
install("tidyverse")
install("tidyfst")
install("igraph")
install("sna")
install("phyloseq")
install("ggalluvial")
install("ggraph")
install("WGCNA")
install("ggnewscale")
install("pulsar")
install("patchwork")
remotes::install_github("taowenmicro/EasyStat")
remotes::install_github("taowenmicro/ggClusterNet")
```


## Example



## 导入R包


```{R}

#--导入所需R包#-------
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(ggClusterNet)

```


## 数据输入

### 内置数据介绍


内置数据为phyloseq格式的数据。该数据来自野生型，突变和过表达三类植物根际，测定的是V5-V7区域的细菌群落；


```{R}
#-----导入数据#-------
data(ps)
ps

```


### 手动数据输入

可以从https://github.com/taowenmicro/R-_function下载数据，构造phylsoeq文件。自己的数据也按照网址示例数据进行准备。虽然phylsoeq对象不易用常规手段处理，但是组学数据由于数据量比较大，数据注释内容多样化，所以往往使用诸如phyloseq这类对象进行处理，并简化数据处理过程。ggClusterNet同样使用了phyloseq对象作为微生物网络的分析。

phyloseq对象构建过程如下，网络分析主要用到otu表格，后续pipeline流程可能用到分组文件metadata，如果按照分类水平山色或者区分模块则需要taxonomy。这几个部分并不是都必须加入phyloseq对象中，可以用到那个加那个。


```{R eval=FALSE, include=FALSE}

library(phyloseq)
library(ggClusterNet)
library(tidyverse)
library(Biostrings)

metadata = read.delim("./metadata.tsv",row.names = 1)
otutab = read.delim("./otutab.txt", row.names=1)
taxonomy = read.table("./taxonomy.txt", row.names=1,header = T)
head(taxonomy)
# tree  = read_tree("./otus.tree")
# rep = readDNAStringSet("./otus.fa")

ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy))#,
              # phy_tree(tree),
              # refseq(rep)
              )

ps
rank.names(ps)

```

### 远程读取数据


但是由于github在国外，所以容易失败

```{R eval=FALSE, include=FALSE}

metadata = read.delim("https://raw.githubusercontent.com/taowenmicro/R-_function/main/metadata.tsv",row.names = 1)
otutab = read.delim("https://raw.githubusercontent.com/taowenmicro/R-_function/main/otutab.txt", row.names=1,header = TRUE)
taxonomy = read.table("https://raw.githubusercontent.com/taowenmicro/R-_function/main/taxonomy.txt", row.names=1,header = TRUE)
head(taxonomy)
head(otutab)
# tree  = read_tree("https://raw.githubusercontent.com/taowenmicro/R-_function/main/otus.tree")
# rep = readDNAStringSet("https://raw.githubusercontent.com/taowenmicro/R-_function/main/otus.fa")

ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy))#,
              # phy_tree(tree),
              # refseq(rep)
              )


```



# 第二章2.1：网络流程新函数network.pip

```{R}

#--更新流程，自由度更高的网络分析流程
library(tidyverse)
library(ggClusterNet)
library(phyloseq)
library(igraph)

#--sparcc方法计算相关矩阵#----
tab.r = network.pip(
  ps = ps,
  N = 200,
  # ra = 0.05,
  big = FALSE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  method = "sparcc",
  label = FALSE,
  lab = "elements",
  group = "Group",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 6
  )


#-提取全部图片的存储对象
plot = tab.r[[1]]
# 提取网络图可视化结果
p0 = plot[[1]]
#zipi
p0.1 = plot[[2]]
#--随机网络幂率分布
p0.2 = plot[[3]]

#--提取相关矩阵,这是一个list存储的相关矩阵
dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab

# 大型相关矩阵跑出来不容易，建议保存，方便各种网络性质的计算
saveRDS(cortab,"cor.matrix.all.group.rds")
cor = readRDS("./cor.matrix.all.group.rds")

```


### 网络显著性分析

```{R}

#--网络显著性比较#-----
dat = module.compare.net.pip(
  ps = NULL,
  corg = cortab,
  degree = TRUE,
  zipi = FALSE,
  r.threshold= 0.8,
  p.threshold=0.05,
  method = "spearman",
  padj = F,
  n = 3)
res = dat[[1]]
head(res)


```

### network.pip 结果如何定制

```{R}
#--网络自定义可视化#-----

dat = tab.r[[2]]
node = dat$net.cor.matrix$node
edge = dat$net.cor.matrix$edge
head(edge)

p <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor
                                  ),
                              data = edge, size = 0.03,alpha = 0.1) +
  geom_point(aes(X1, X2,
                 fill = Phylum,
                 size = igraph.degree),
             pch = 21, data = node,color = "gray40") +
  facet_wrap(.~ label,scales="free_y",nrow = 1) +
  # geom_text_repel(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
  # geom_text(aes(X1, X2,label = elements),pch = 21, data = nodeG) +
  scale_colour_manual(values = c("#6D98B5","#D48852")) +
  # scale_fill_hue()+
  scale_size(range = c(0.8, 5)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()
  ) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p

```

### zipi结果定制

```{R}

#--zipi可视化-定制#-----
dat.z = dat$zipi.data
head(dat.z)
x1<- c(0, 0.62,0,0.62)
x2<- c( 0.62,1,0.62,1)
y1<- c(-Inf,2.5,2.5,-Inf)
y2 <- c(2.5,Inf,Inf,2.5)
lab <- c("peripheral",'Network hubs','Module hubs','Connectors')
roles.colors <- c("#E6E6FA","#DCDCDC","#F5FFFA", "#FAEBD7")
tab = data.frame(x1 = x1,y1 = y1,x2 = x2,y2 = y2,lab = lab)
tem = dat.z$group %>% unique() %>% length()
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
  geom_point(data=dat.z,aes(x=p, y=z,color=module)) + theme_bw()+
  guides(color= F) +
  ggrepel::geom_text_repel(data = dat.z,
                           aes(x = p, y = z,
                               color = module,label=label),size=4)+
  # facet_wrap(.~group) +
  facet_grid(.~ group, scale='free') +
  theme(strip.background = element_rect(fill = "white"))+
  xlab("Participation Coefficient")+ylab(" Within-module connectivity z-score")
p
```

### 网络与随机网络比对

```{R}

# --随机网络，幂率分布#-------
dat.r = dat$random.net.data

p3 <- ggplot(dat.r) +
  geom_point(aes(x = ID,y = network,
                 group =group,fill = group),pch = 21,size = 2) +
  geom_smooth(aes(x = ID,y = network,group =group,color = group))+
  facet_grid(.~g,scales = "free") +
  theme_bw() + theme(
    plot.margin=unit(c(0,0,0,0), "cm")
  )
p3
```

# Additional Information
## ggClusterNet

Microbial ecological network visualization clustering

![](http://www.imeta.science/iMeta/Papers/8GraphicAbstract/imt2.32.jpg)

Citation: Tao Wen, Penghao Xie, Shengdie Yang, Guoqing Niu, Xiaoyu Liu, Zhexu Ding, Chao Xue, Yong-Xin Liu, Qirong Shen, Jun Yuan. 2022. ggClusterNet: An R package for microbiome network analysis and modularity-based multiple network layouts. iMeta 1: e32. https://doi.org/10.1002/imt2.32


