#' Microbial matel text and env correcton with the plot to showing,from science peper.
#'
#' @param env.dat env table
#' @param tabOTU OTU table conbined with list
#' @param distance The defult, distance = "bray", method of microbiota distance
#' @param p.threshold The defult, p.threshold=0.05, it represents significance threshold below 0.05.
#' @param method method for matel test
#' @param x x axis label
#' @param y y axis label
#' @param diag diag label
#' @param fill fill coulor of node
#' @param sig if show only sig value
#' @param siglabel siglabel
#' @param siglabel R2 value
#' @param numpoint point type 22
#' @param numpoint2 point type 21
#' @param numsymbol point adding  e.g. 27
#' @param lacx location of x axis
#' @param lacy location of y axis
#' @param range curve size
#' @param p.thur sig value
#' @param onlysig T or F
#' @examples
#' data(ps)
#' library(phyloseq)
#' library(tidyverse)
#' ps1 = filter_taxa(ps, function(x) sum(x ) > 200 , TRUE);ps1
#' otu = as.data.frame(t(vegan_otu(ps1)))
#' mapping = as.data.frame( sample_data(ps1))
#' tabOTU1 = list(otu1 = otu,otu2 = otu,otu3 = otu)
#' data(env1)
#' MatCorPlot(
#'   env.dat = env1,
#'  tabOTU = tabOTU1,
#'   distance = "bray",
#'   method = "metal",
#'   x = F,
#'   y = F,
#'   diag = T,
#'   sig = TRUE,
#'   siglabel = FALSE,
#'   shownum = TRUE,
#'   numpoint = NULL,
#'   numsymbol = 27,
#'   lacx = "left",
#'   lacy = "bottom",
#'   range = 0.5,
#'   p.thur = 0.3,
#'   onlysig = F
#' )
#' @return ggplot object
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export

#---tagether three sub function
MatCorPlot <- function(
  env.dat,
  tabOTU,
  distance = "bray",
  method = "metal",
  method.cor = "spearman",
  cor.p = 0.05,
  x = TRUE,
  y = TRUE,
  diag = T,
  sig = TRUE,
  siglabel = FALSE,
  shownum = TRUE,
  numpoint = 22,
  numpoint2 = 21,
  numsymbol = NULL,
  curvature = 0.2,
  lacx = "left",
  lacy = "bottom",
  range = 0.5,# set of line size
  p.thur = 0.3,# sig value
  onlysig = TRUE# if show sig connect
){

  #--- mantel test
  rep = MetalTast (env.dat = env.dat, tabOTU = tabOTU,distance = distance,method = method)
  repR = rep[c(-seq(from=1,to=dim(rep)[2],by=2)[-1])]
  repP = rep[seq(from=1,to=dim(rep)[2],by=2)]


  result <- Miccorplot(data = env.dat,
                       method.cor = method.cor,
                       cor.p = cor.p,
                       x = x,
                       y = y,
                       diag = diag,
                       lacx = lacx,
                       lacy = lacy,
                       sig = sig,
                       siglabel = siglabel,
                       shownum = shownum,
                       numpoint = numpoint,
                       numsymbol = numsymbol)


  cor_link(
    data = result[[2]],# env table
    p = result[[1]], # ggplot object
    envdata = repR,# R value table
    Ptab = repP,#p value table
    numpoint2 = numpoint2,
    curvature = curvature,
    range = range,# line size
    lacx = lacx,#  left or right
    lacy = lacy,# top or bottom
    p.thur = p.thur,#sig value
    onlysig = onlysig# if show sig connect
  )
}


#---function for cor matrix of helf
Miccorplot <- function(data,
                       method.cor="spearman",
                       cor.p = 0.05,
                       x = TRUE,
                       y = TRUE,
                       diag = TRUE,
                       lacx = "left",
                       lacy = "bottom",
                       sig = TRUE,
                       siglabel = FALSE,
                       shownum = TRUE,
                       numpoint = 22,
                       numsymbol = NULL,
                       dig.vjust = -0.5,
                       dig.hjust = -0.5
){

  lab.1 = 0
  lab.2 = 0
  #-culculate cor
  occor = psych::corr.test(data,use="pairwise",method=method.cor,adjust="fdr",alpha=cor.p)
  occor.r = occor$r
  occor.p = occor$p
  if (sig == TRUE) {
    occor.r[occor.p > cor.p] = 0
  }
  #--pick out lower matrix
  occor.r2 = occor.r[lower.tri(occor.r, diag = TRUE)]
  #--cul x,y
  n = dim(occor.r)[1]
  axis = 1:n
  axis.2 = axis[order(axis, decreasing=TRUE)]
  #--lab#---
  colname = colnames(occor.r)
  rowname = row.names(occor.r)

  #--left#----
  if (lacx == "left") {
    axis.x = matrix(1:n,n,n,byrow = TRUE)[lower.tri(matrix(axis,n,n,byrow = TRUE),diag = TRUE)]
    lab.1 = lab.1
  }
  # right
  if (lacx == "right") {
    axis.x = matrix(n:1,n,n,byrow = TRUE)[lower.tri(matrix(axis,n,n,byrow = TRUE),diag = TRUE)]
    lab.1 = n+1
  }

  #--bottom#----
  if (lacy == "bottom") {
    axis.y = matrix(n:1,n,n,byrow = FALSE)[lower.tri(matrix(n:1,n,n,byrow = FALSE),diag = TRUE)]
    lab.2 = lab.2
  }
  #--top
  if (lacy == "top") {
    axis.y = matrix(1:n,n,n,byrow = FALSE)[lower.tri(matrix(1:n,n,n,byrow = FALSE),diag = TRUE)]
    lab.2 = n+1
  }
  # plot data
  df0 = data.frame(x = axis.x,y = axis.y,r = occor.r2)
  labdat <- data.frame(x = axis,z = axis.2,labx = colname,laby = colname[axis.2] )
  # base plot
  p <- ggplot() +
    geom_tile(aes(x = x, y = y),data = df0,fill=NA,color='gray',size=0.5) +
    scale_size(range = c(1, 8))+
    scale_fill_distiller(palette="RdYlBu") + theme_void()
  #--if sig only showing
  if (sig == TRUE) {
    df = df0 %>% filter(r != 0)
    addat = df
  }  else if(sig == FALSE){
    addat <- df0
  }

  # p <- p + geom_point(aes(x = x, y = y,size=r,fill=r),df0, shape=22,color='white')
  if (!is.null(numpoint) ) {
    p <- p + geom_point(aes(x = x, y = y,size=r,fill=r),addat, shape=numpoint,color='black')
  }else if (!is.null(numsymbol) & is.null(numpoint)) {
    p <- p  + ggsymbol:: geom_symbol(data = addat,aes(x = x, y = y,size=r,fill = r),symbol = numsymbol)
  }


  # ig dig only  adding *
  if (siglabel == TRUE) {
    addat$r2 = "*"
    p <- p + geom_text(aes(x = x, y = y,label = r2),addat)

  }

  if (shownum == TRUE) {
    p <- p + geom_text(aes(x = x, y = y,label = round(r,2)),addat)
  }

  # adding label#------
  if (x == TRUE) {
    p <- p + geom_text(aes(x = x, y = lab.2,label= labx),labdat)
  }

  if (y == TRUE) {
    p <- p + geom_text(aes(x = lab.1, y = x,label= laby),labdat)
  }


  if (lab.1 == 0 ) {
    if (diag == TRUE) {
      p <- p +  geom_text(aes(x = x + 1, y = abs(lab.2-z),label=labx),labdat,angle = 30,vjust =  dig.vjust,hjust = dig.vjust)
    }
    lindat = data.frame(x = labdat$x + 1,y = abs(lab.2-labdat$z) ,label = labdat$labx)
  }

  if (lab.1 == (n+1) & lab.2 ==0) {
    if (diag == TRUE) {
      p <- p + geom_text(aes(x = x - 1, y = abs(lab.1-z),label=labx),labdat,angle = -30,vjust =  -0.5)
    }
    lindat = data.frame(x = labdat$x - 1,y = abs(lab.1-labdat$z),label = labdat$labx)
  } else if (lab.1 == (n+1) & lab.2 == (n+1)) {
    if (diag == TRUE) {
      p <- p + geom_text(aes(x = x - 1, y = abs(z),label=labx),labdat,angle = 30,vjust =  1)
    }

    lindat = data.frame(x = labdat$x - 1,y = abs(labdat$z),label = labdat$labx )
  }


  return(list(p,lindat))
}



cor_link <- function(data,
                     p,
                     envdata = report,
                     Ptab,
                     range = 1,
                     lacx = "left",
                     lacy = "bottom",
                     numpoint2 = 21,
                     curvature = 0.2,
                     p.thur = 0.3,
                     onlysig = TRUE
){

  colnames(envdata)[1] = c("label")
  data3<- data %>% inner_join(envdata)

  colnames(Ptab)[1] = "label"
  data3<- data3 %>% inner_join(Ptab)


  for (i in 1:length(colnames(envdata)[-1])) {
    group <- c()
    group[data3[colnames(envdata)[-1][i]] > 0] <- "pos"
    group[data3[colnames(envdata)[-1][i]] < 0] <- "neg"
    group[data3[colnames(envdata)[-1][i]] == 0] <- "nose"

    group <- data.frame(label= data3$label ,group)
    colnames(group)[2] <- paste("group",colnames(envdata)[-1][i],sep = "")
    data3 <- data3 %>% inner_join(group)
  }

  #-p matrix
  for (i in 1:length(colnames(Ptab)[-1])) {
    group <- c()
    group[data3[colnames(Ptab)[-1][i]] > p.thur] <- "nose"
    group[data3[colnames(Ptab)[-1][i]] < p.thur] <- "sig"

    group <- data.frame(label= data3$label ,group)
    colnames(group)[2] <- paste("group",colnames(Ptab)[-1][i],sep = "")
    data3 <- data3 %>% inner_join(group)
  }

  # filter R2 for curve
  if (onlysig == TRUE) {
    for (i in 1:length(colnames(Ptab)[-1])) {
      data3[colnames(envdata)[-1][i]][data3[ paste("group",colnames(Ptab)[-1][i],sep = "")] == "nose"] = NA

    }
  }




  if (lacx == "left" &lacy == "bottom") {
    linx = c(max(data3$x)*0.7,max(data3$x),max(data3$x)*1.3)
  } else if(lacx == "left" &lacy == "top") {
    linx = c(max(data3$x)*0.35,max(data3$x)*0.75,max(data3$x)*1.05)
  }else if(lacx == "right" &lacy == "top") {
    linx = c(max(data3$x)*0.35,-max(data3$x)*0.35,-max(data3$x)*0.75)
  }else if(lacx == "right" &lacy == "bottom") {
    linx = c(-max(data3$x)*0.75,-max(data3$x)*0.35,0)
  }

  if (lacx == "left" &lacy == "bottom") {
    liny = c(max(data3$y)*1.3,max(data3$y) ,max(data3$y)*0.7)

  } else if (lacx == "left" &lacy == "top") {
    liny = c(max(data3$x)*0.05,max(data3$x)*0.25,max(data3$y)*0.7)
  }else if (lacx == "right" &lacy == "top") {
    liny = c(max(data3$x)*0.05,max(data3$x)*0.25,max(data3$y)*0.7)
  }else if (lacx == "right" &lacy == "bottom") {
    liny = c(max(data3$x)*0.25,max(data3$x)*0.75,max(data3$y)*1.2)
  }


  n = length(colnames(envdata)[-1])

  if (n == 1 & lacx == "left") {

    p = p +
      geom_curve(data = data3,aes_string(x = max(data3$x), y = max(data3$x)*0.5, xend = "x", yend = "y",
                                         group = paste("group",colnames(envdata)[-1][1],sep = ""),
                                         color = paste("group",colnames(envdata)[-1][1],sep = "")),
                 curvature = curvature,size = data3[,colnames(envdata)[-1][1]]^2*200 *range) +
      geom_point(data = data3,aes(x = x, y = y),pch = numpoint2,size =4,color = "black",fill = "#FFF5EB")+
      geom_point(data = data3,aes(x = max(data3$x), y = max(data3$x)*0.5),pch = numpoint2,size = 6,color = "black",fill = "#FEE6CE") +
      geom_text(data = data3,aes(x = max(data3$x), y =max(data3$x)*0.5,label = colnames(envdata)[-1][1] ),vjust = -1,hjust = -1)


  } else if(n == 1 & lacx == "right") {
    p = p +
      geom_curve(data = data3,aes_string(x = 0, y = max(data3$x)*0.5, xend = "x", yend = "y",
                                         group = paste("group",colnames(envdata)[-1][1],sep = ""),
                                         color = paste("group",colnames(envdata)[-1][1],sep = "")),
                 curvature = curvature,size = data3[,colnames(envdata)[-1][1]]^2*200 *range) +
      geom_point(data = data3,aes(x = x, y = y),pch = numpoint2,size =4,color = "black",fill = "#FFF5EB")+
      geom_point(data = data3,aes(x = 0, y = max(data3$x)*0.5),pch = numpoint2,size = 6,color = "black",fill = "#FEE6CE") +
      geom_text(data = data3,aes(x = 0, y =max(data3$x)*0.5,label = colnames(envdata)[-1][1] ),vjust = -1,hjust = -1)

  }



  if (n == 2) {
    p = p +
      geom_curve(data = data3,aes_string(x =linx[1], y = liny[1], xend = "x", yend = "y",group = paste("group",colnames(envdata)[-1][1],sep = ""),
                                         color = paste("group",colnames(envdata)[-1][1],sep = "")), curvature = curvature,size = data3[,colnames(envdata)[-1][1]]^2*200*range) +
      geom_curve(data = data3,aes_string(x =linx[3], y = liny[3], xend = "x", yend = "y",group = paste("group",colnames(envdata)[-1][2],sep = ""),
                                         color = paste("group",colnames(envdata)[-1][2],sep = "")), curvature = curvature,size = data3[,colnames(envdata)[-1][2]]^2*200*range) +
      geom_point(data = data3,aes(x = x, y = y),pch = numpoint2,size =4,color = "black",fill = "#FFF5EB")+
      geom_point(data = data3,aes(x =linx[1], y = liny[1]),pch = numpoint2,size = 6,color = "black",fill = "#FEE6CE") +
      geom_point(data = data3,aes(x =linx[3], y = liny[3]),pch = numpoint2,size = 6,color = "black",fill = "#FEE6CE") +
      geom_text(data = data3,aes(x =linx[1], y = liny[1],label = colnames(envdata)[-1][1] ),vjust = -1,hjust = -1) +
      geom_text(data = data3,aes(x =linx[3], y = liny[3],label = colnames(envdata)[-1][2] ),vjust = -1,hjust = -1)
  }



  if (n == 3) {
    p = p +
      geom_curve(data = data3,aes_string(x =linx[1], y = liny[1], xend = "x", yend = "y",group = paste("group",colnames(envdata)[-1][1],sep = ""),
                                         color = paste("group",colnames(envdata)[-1][1],sep = "")), curvature = curvature,size = data3[,colnames(envdata)[-1][1]]^2*200*range) +
      geom_curve(data = data3,aes_string(x =linx[2], y = liny[2], xend = "x", yend = "y",group = paste("group",colnames(envdata)[-1][2],sep = ""),
                                         color = paste("group",colnames(envdata)[-1][2],sep = "")), curvature = curvature,size = data3[,colnames(envdata)[-1][2]]^2*200*range) +
      geom_curve(data = data3,aes_string(x =linx[3], y = liny[3], xend = "x", yend = "y",group = paste("group",colnames(envdata)[-1][3],sep = ""),
                                         color = paste("group",colnames(envdata)[-1][3],sep = "")), curvature = curvature,size = data3[,colnames(envdata)[-1][3]]^2*200*range) +
      geom_point(data = data3,aes(x = x, y = y),pch = numpoint2,size =4,color = "black",fill = "#FFF5EB") +
      geom_point(data = data3,aes(x =linx[1], y = liny[1]),pch = numpoint2,size = 6,color = "black",fill = "#FEE6CE") +
      geom_point(data = data3,aes(x =linx[2], y = liny[2]),pch = numpoint2,size = 6,color = "black",fill = "#FEE6CE") +
      geom_point(data = data3,aes(x =linx[3], y = liny[3]),pch = numpoint2,size = 6,color = "black",fill = "#FEE6CE") +
      geom_text(data = data3,aes(x =linx[1], y = liny[1],label = colnames(envdata)[-1][1] ),vjust = -1,hjust = -1) +
      geom_text(data = data3,aes(x =linx[2], y = liny[2],label = colnames(envdata)[-1][2] ),vjust = -1,hjust = -1) +
      geom_text(data = data3,aes(x =linx[3], y = liny[3],label = colnames(envdata)[-1][3] ),vjust = -1,hjust = -1)

  }

  # colnames(Ptab)[-1][1]


  p
}

# env.dat = env
# tabOTU = tabOTU1
# MetalTast(env.dat = env,tabOTU = tabOTU1,distance = "jcd",method = "metal")


MetalTast <- function(env.dat,
                      tabOTU,
                      distance = "bray",
                      method = "metal"
    ){

  for (j in 1:length(names(tabOTU))) {
    report =c()
    otu = tabOTU[[j]]
    name = names(tabOTU)[j]
    sle_env  = vegan::decostand(env.dat, method="standardize", MARGIN=2)#
    ##env distance
    env.std = vegan::decostand(sle_env, method = "standardize", MARGIN=2) #按列标准化到均值等于0，标准偏差等于1.

    ##OTU distance
    BC.beta = vegan::vegdist(t(otu), method="bray")
    JC.beta = vegan::vegdist(t(otu), method="jaccard",binary=T)
    if (method == "metal") {
      for(i in 1:ncol(env.std)){
        envdis =vegan::vegdist(env.std[i],method = "euclidean", na.rm=T)
        mantel.BC = vegan::mantel(envdis, BC.beta, na.rm=T)
        mantel.JC = vegan::mantel(envdis, JC.beta, na.rm=T)
        report = rbind(report,c(colnames(env.dat)[i], mantel.BC$statistic, mantel.BC$signif, mantel.JC$statistic, mantel.JC$signif))
      }

    }

    if (method == "Part.Mantel") {
      for(i in 1:ncol(env.std)){   ##
        envdis = dist(env.std[,i])
        envdis2 = dist(env.std[,-i])
        pmantel.BC = mantel.partial(BC.beta, envdis, envdis2, na.rm=T)
        pmantel.JC = mantel.partial(JC.beta, envdis, envdis2, na.rm=T)
        report = rbind(report,c(colnames(env.dat)[i], pmantel.BC$statistic, pmantel.BC$signif, pmantel.JC$statistic, pmantel.JC$signif))
      }
    }

    colnames(report)<- c("Envs",paste(rep(c("r","p"),2),rep(c("BC","JC"),each=2),sep = "."))
    report = as.matrix(report)
    report = as.data.frame(report)
    report[,2:dim(report)[2]]<-lapply(report[,2:dim(report)[2]],as.character)
    report[,2:dim(report)[2]]<-lapply(report[,2:dim(report)[2]],as.numeric)



    #首先我们使用BARY距离计算的mantel结果进行计算



    if (distance == "bray" ) {
      report0 = report[1:3]
      colnames(report0)[-1] = paste(name,colnames(report0)[-1],sep = "")
    } else if (distance == "jcd"){
      report0 = report[c(1,c(4:5))]
      colnames(report0)[-1] = paste(name,colnames(report0)[-1],sep = "")
    }

    if (j == 1) {
      report1 = report0
    }

    if (j != 1) {
      report1 = report1 %>% inner_join(report0)
    }

  }

  return(report1)
}

