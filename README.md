

#  ggClusterNet 2.0: an R package for microbial co-occurrence networks and associated indicator correlation patterns

Since the last version release in 2022, ggClusterNet has emerged as a critical resource for microbiome research, enabling microbial co-occurrence network analysis and visualization in over 200 studies (Google Scholar citations). To address emerging challenges in microbiome studies, including multi-factor experimental designs, multi-treatment, and multi-omics data, we present a comprehensive upgrade with the following four components: 1) We recommended and designed a microbial co-occurrence network analysis pipeline incorporating network computation and visualization (Pearson/Spearman/SparCC correlations), topological characterization of network and node properties, multi-network structure comparison and statistical testing, exploration of network stability (robustness), and identification and analysis of network modules; 2) Developed microbial network mining functions for multi-factor, multi-treatment, and spatiotemporal-scale analysis, such as Facet. Network(), module.compare.m.ts(), Robustness.Random.removal.ts(), etc.; 3) Developed functions for microbial and multi-factor interaction analysis, along with versatile visualization layout algorithms, such as MatCorPlot2(), Miccorplot3(), cor_link3(), matcorplotj(), and two.cor(); 4) Developed functions for cross-domain and multi-omics integrated network analysis, including corBionetwork.st(), and developed a comprehensive suite of visualization layout algorithms specifically designed for exploring complex relationships in these networks, such as model_maptree2(), model_Gephi.3(), cir.squ(), and cir.maptree2(). Collectively, the latest updates to ggClusterNet 2.0 empower researchers to explore complex network interactions with enhanced capabilities, offering a robust, efficient, user-friendly, reproducible, and visually versatile tool for microbial co-occurrence networks and associated indicator correlation patterns. The ggClusterNet 2.0 R package is open-source and freely accessible on GitHub (https://github.com/taowenmicro/ggClusterNet).



## Examples of visualizations.

![](https://github.com/taowenmicro/Rcoding/blob/main/ggClusternet2.0.plot/Fig3.jpg?raw=true)


![](https://github.com/taowenmicro/Rcoding/blob/main/ggClusternet2.0.plot/Fig4.jpg?raw=true)


## Main features:
- 1)	The ggClusterNet 2 introduces a comprehensive microbial co-occurrence network analysis pipeline.
- 2)	Enhances the network analysis workflow tailored for complex experimental designs and diverse data types.
- 3)	Enhances visualization capabilities for exploring microbiomes and their correlated environmental or host-associated indicators.
- 4)	Introduces a variety of visualization layout algorithms suitable for cross-domain and multi-omics interaction networks.



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

## input data

###  data

phyloseq

```{R}

data(ps)
ps

```

# network.pip

```{R}


library(tidyverse)
library(ggClusterNet)
library(phyloseq)
library(igraph)

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


plot = tab.r[[1]]
p0 = plot[[1]]
p0
p0.1 = plot[[2]]
p0.2 = plot[[3]]

dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab

saveRDS(cortab,"cor.matrix.all.group.rds")
cor = readRDS("./cor.matrix.all.group.rds")

```


##  Reference

If used this script, please cited:

Tao Wen, Penghao Xie, Shengdie Yang, Guoqing Niu, Xiaoyu Liu, Zhexu Ding, Chao Xue, Yong-Xin Liu, Qirong Shen, Jun Yuan. 2022. ggClusterNet: An R package for microbiome network analysis and modularity-based multiple network layouts. iMeta 1: e32. https://doi.org/10.1002/imt2.32



