


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
#' MatCorPlot2(
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
MatCorPlot2 <- function(
    env.dat,
    tabOTU,
    distance = "bray",
    method = "metal",
    method.cor = "spearman",
    cor.p = 0.05,
    x = TRUE,
    y = TRUE,
    zoom = 2,
    corva = 0.2,
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
  rep = MetalTast(env.dat = env.dat, tabOTU = tabOTU,distance = distance,method = method)
  repR = rep[c(-seq(from=1,to=dim(rep)[2],by=2)[-1])]
  repP = rep[seq(from=1,to=dim(rep)[2],by=2)]


  result <- Miccorplot2(data = env.dat,
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



  cor_link2(
    data = result[[2]],# env table
    p = result[[1]], # ggplot object
    envdata = repR,# R value table
    Ptab = repP,#p value table
    zoom = zoom,
    corva = corva,
    numpoint2 = numpoint2,
    curvature = curvature,
    range = range,# line size
    lacx = lacx,#  left or right
    lacy = lacy,# top or bottom
    p.thur = p.thur,#sig value
    onlysig = onlysig # if show sig connect
  )
}


#---function for cor matrix of helf
Miccorplot2 <- function(data,
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
  occor = psych::corr.test(data,
                           use="pairwise",
                           method=method.cor,
                           adjust="fdr",alpha=cor.p)
  occor.r = occor$r
  occor.p = occor$p
  if (sig == TRUE) {
    occor.r[occor.p > cor.p] = 0
  }
  #--pick out lower matrix
  occor.r2 = occor.r[lower.tri(occor.r, diag = FALSE)]
  #--cul x,y
  n = dim(occor.r)[1]
  axis = 1:n
  axis.2 = axis[order(axis, decreasing=TRUE)]
  #--lab#---
  colname = colnames(occor.r)
  rowname = row.names(occor.r)

  #--left#----
  n = n -1
  if (lacx == "left") {
    axis.x = matrix(1:n,n,n,byrow = TRUE)[lower.tri(matrix(axis,n,n,byrow = TRUE),
                                                    diag = TRUE)]
    lab.1 = lab.1
  }
  # right
  if (lacx == "right") {
    axis.x = matrix(n:1,n,n,byrow = TRUE)[lower.tri(matrix(axis,n,n,byrow = TRUE),
                                                    diag = TRUE)]
    lab.1 = n+1
  }

  #--bottom#----
  if (lacy == "bottom") {
    axis.y = matrix(n:1,n,n,byrow = FALSE)[lower.tri(matrix(n:1,n,n,byrow = FALSE),
                                                     diag = TRUE)]
    lab.2 = lab.2
  }
  #--top
  if (lacy == "top") {
    axis.y = matrix(1:n,n,n,byrow = FALSE)[lower.tri(matrix(1:n,n,n,byrow = FALSE),diag = TRUE)]
    lab.2 = n+1
  }


  # plot data
  df0 = data.frame(x = axis.x,y = axis.y,r = occor.r2)
  labdat <- data.frame(x = axis,z = axis.2,
                       labx = colname[length(colname):1],
                       laby = colname[axis.2])
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
  head(addat)
  addat$r2 = abs(addat$r)
  colnames(addat)[4] = paste0(stringr::str_to_title(method.cor),".r")
  colnames(addat)[3] = paste0(stringr::str_to_title(method.cor),"_r")
  # p <- p + geom_point(aes(x = x, y = y,size=r,fill=r),df0, shape=22,color='white')
  if (!is.null(numpoint)) {
    p <- p + geom_point(aes(x = x, y = y,
                            size=!!sym(colnames(addat)[4]),fill=!!sym(colnames(addat)[3])),addat, shape=numpoint,color='black')+
      scale_fill_gradient2(low = "#008B45FF",
                           high = "#EE0000FF") +
      theme(plot.margin=unit(rep(2,4),'cm'))

  }else if (!is.null(numsymbol) & is.null(numpoint)) {
    p <- p  +
      ggsymbol:: geom_symbol(data = addat,aes(x = x, y = y,size=!!sym(colnames(addat)[4]),
                                              fill = !!sym(colnames(addat)[3])),
                             symbol = numsymbol)+
      theme(plot.margin=unit(rep(2,4),'cm'))
  }


  # ig dig only  adding *
  if (siglabel == TRUE) {
    addat$r2 = "*"
    p <- p + geom_text(aes(x = x, y = y,label = r2),addat)

  }

  if (shownum == TRUE) {
    p <- p + geom_text(aes(x = x, y = y,label = round(!!sym(colnames(addat)[3]),2)),addat)
  }

  # adding x label#------
  if (x == TRUE&lacx == "right") {
    p <- p + geom_text(aes(x = x -1, y = lab.2,label= labx),labdat[-1,])+
      theme(plot.margin=unit(rep(2,4),'cm'))
  } else if(x == TRUE&lacx == "left"){
    labdat$x2 = (length(labdat$x) - labdat$x)
    p <- p + geom_text(aes(x = x2 +1, y = lab.2,label= labx),labdat[-1,])+
      theme(plot.margin=unit(rep(2,4),'cm'))
  }
  # add y label#-------
  if (y == TRUE & lacx == "right" & lacy == "bottom") {
    p <- p + geom_text(aes(x = lab.1, y = x,label= laby),labdat[-nrow(labdat),])+
      theme(plot.margin=unit(rep(2,4),'cm'))
  }else if(y == TRUE&lacx == "left"&lacy == "bottom"){
    labdat$x2 = length(labdat$x) - labdat$x
    p <- p + geom_text(aes(x = lab.1, y = x2,label= laby),labdat[-nrow(labdat),])+
      theme(plot.margin=unit(rep(2,4),'cm'))
  } else if (y == TRUE & lacx == "right"&lacy == "top"){
    labdat$x2 = length(labdat$x) - labdat$x
    p <- p + geom_text(aes(x = lab.1, y = x2,label= laby),labdat[-nrow(labdat),])+
      theme(plot.margin=unit(rep(2,4),'cm'))
  }else if(y == TRUE&lacx == "left"&lacy == "top"){
    p <- p + geom_text(aes(x = lab.1, y = x,label= laby),labdat[-nrow(labdat),])+
      theme(plot.margin=unit(rep(2,4),'cm'))
  }


  if (lab.1 == 0 ) {
    if (diag == TRUE&lacy == "bottom"&lacx == "right") {
      labdat$x = length(labdat$x) - labdat$x
      labdat$z = length(labdat$x) - labdat$z
      p <- p +  geom_text(aes(x = x+1 , y = z +1,label=labx),labdat,
                          vjust =  1,hjust = 0.5)+
        theme(plot.margin=unit(rep(2,4),'cm'))
      lindat = data.frame(x = labdat$x + 1,y = labdat$z +1 ,label = labdat$labx)
    } else if(diag == TRUE&lacy == "top"&lacx == "right"){
      # labdat$x = length(labdat$x) - labdat$x
      # labdat$z = length(labdat$x) - labdat$z
      p <- p +  geom_text(aes(x = x+1 , y = z -1,label=labx),labdat,
                          vjust =  1,hjust = 0.5)+
        theme(plot.margin=unit(rep(2,4),'cm'))
      lindat = data.frame(x = labdat$x + 1,y = labdat$z - 1 ,label = labdat$labx)
    }  else if(diag == TRUE&lacy == "bottom"&lacx == "left") {
      labdat$x = length(labdat$x) - labdat$x
      labdat$z = length(labdat$x) - labdat$z
      p <- p +  geom_text(aes(x = x +1, y = z+1 ,label=labx),labdat,
                          vjust =  1,hjust = 0.5)+
        theme(plot.margin=unit(rep(2,4),'cm'))
      lindat = data.frame(x = labdat$x + 1,y = labdat$z +1 ,label = labdat$labx)
    } else if(diag == TRUE&lacy == "top"&lacx == "left"){
      labdat$x = length(labdat$x) - labdat$x
      # labdat$z = length(labdat$x) - labdat$z
      p <- p +  geom_text(aes(x = x+1 , y = z -1,label=labx),labdat,
                          vjust =  1,hjust = 0.5)+
        theme(plot.margin=unit(rep(2,4),'cm'))
      lindat = data.frame(x = labdat$x + 1,y = labdat$z - 1 ,label = labdat$labx)
    }




  }


  if (lab.1 == (n+1) & lab.2 ==0) {
    if (diag == TRUE) {
      p <- p + geom_text(aes(x = x-1 , y = abs(lab.1-z+1),label=labx),labdat,vjust =  1,hjust = 0.5)+
        theme(plot.margin=unit(rep(2,4),'cm'))
    }

    lindat = data.frame(x = labdat$x - 1,y = lab.1- labdat$z +1,label = labdat$labx)

  } else if (lab.1 == (n+1) & lab.2 == (n+1)) {
    if (diag == TRUE) {
      labdat$x2 = length(labdat$x) - labdat$x
      # labdat$z = length(labdat$x) - labdat$z
      p <- p + geom_text(aes(x = x -1  , y = z-1,label=labx),labdat,vjust =  1,hjust = 0.5)+
        theme(plot.margin=unit(rep(2,4),'cm'))
    }

    lindat = data.frame(x = labdat$x - 1,y = labdat$z -1 ,label = labdat$labx )
  }


  return(list(p,lindat))
}



cor_link2 <- function(data,
                     p,
                     envdata = report,
                     Ptab,
                     zoom = 2,
                     corva = 0.2,
                     range = 1,
                     lacx = "left",
                     lacy = "bottom",
                     numpoint2 = 21,
                     curvature = 0.2,
                     p.thur = 0.3,
                     onlysig = TRUE
){


  colnames(envdata)[1] = c("label")
  dat1 = envdata %>%gather(key = "groupbc", value = "R",-label)

  data31<- data %>% full_join(dat1,by = "label")

  colnames(Ptab)[1] = "label"
  dat2 = Ptab %>%gather(key = "group.p.bc", value = "p",-label)
  head(dat2)
  data32 <- data %>% left_join(dat2,by = "label")
  data3 = cbind(data31,data32[,-(1:3)])
  head(data3)

    group <- c()
    group[data3$R > 0] <- "pos"
    group[data3$R < 0] <- "neg"
    group[data3$R == 0] <- "nose"
    data3$groupr = group


  #-p matrix
    group <- c()
    group[data3$p > p.thur] <- "nose"
    group[data3$p < p.thur] <- "sig"

    data3$groupp = group
  # head(data3)
    data3$group = gsub("p.BC","",data3$group.p.bc)
  x = length(unique(data3$group))
  idtab = data3
  # filter R2 for curve
  if (onlysig == TRUE) {
    data3 = data3 %>% filter(p < p.thur )
  }

  # ggplot() + geom_point(aes(x,y),data = data) +
  #   geom_point(aes(x -1,y +1),data = data)
  #
  # ggplot() + geom_point(aes(x,y),data = data) +
  #   geom_point(aes(xend,yend),data = topdat)
  n = x +1
  if (lacx == "left"& lacy == "bottom") {
    dat = data.frame(x = data$x + zoom,y = data$y + zoom)
    seqnum = (max(dat$x) - min(dat$x))/n
    tem = seq(from=min(dat$y), to=max(dat$y),by=seqnum)
    tem = tem[length(tem):1]

    topdat = data.frame(x=seq(from=min(dat$x), to=max(dat$x),by=seqnum),
                        y = tem) %>% slice(-1)
    topdat = topdat[-nrow(topdat),]
    colnames(topdat) = paste0(colnames(topdat),"end")
    topdat$group = unique(data3$group)

  } else if (lacx == "right"& lacy == "bottom"){
    dat = data.frame(x = data$x - zoom,y = data$y + zoom)
    seqnum = (max(dat$x) - min(dat$x))/n
    topdat = data.frame(x=seq(from=min(dat$x), to=max(dat$x),by=seqnum),
                        y = seq(from=min(dat$y), to=max(dat$y),by=seqnum)) %>% slice(-1)
    topdat = topdat[-nrow(topdat),]
    colnames(topdat) = paste0(colnames(topdat),"end")
    topdat$group = unique(data3$group)
  } else if(lacx == "right"&lacy == "top"){
    dat = data.frame(x = data$x - zoom,y = data$y - zoom)
    seqnum = (max(dat$x) - min(dat$x))/n
    tem = seq(from=min(dat$y), to=max(dat$y),by=seqnum)
    tem = tem[length(tem):1]

    topdat = data.frame(x=seq(from=min(dat$x), to=max(dat$x),by=seqnum),
                        y = tem) %>% slice(-1)
    topdat = topdat[-nrow(topdat),]
    colnames(topdat) = paste0(colnames(topdat),"end")
    topdat$group = unique(data3$group)

  } else if( lacx == "left"&lacy == "top"){
    dat = data.frame(x = data$x + zoom,y = data$y - zoom)
    seqnum = (max(dat$x) - min(dat$x))/n
    tem = seq(from=min(dat$y), to=max(dat$y),by=seqnum)
    # tem = tem[length(tem):1]

    topdat = data.frame(x=seq(from=min(dat$x), to=max(dat$x),by=seqnum),
                        y = tem) %>% slice(-1)
    topdat = topdat[-nrow(topdat),]
    colnames(topdat) = paste0(colnames(topdat),"end")
    topdat$group = unique(data3$group)

  }





  # head(data3)


  tab2 = data3 %>% left_join(topdat,by = "group")
  # head(tab2)

  # corva = -0.2
  tem = c(rep(corva,floor(nrow(data)/2)),
  rep(-corva,ceiling(nrow(data)/2)))

  # range = 0.02
  tem2 = tab2$R^2*200 *range
  tab2$size = tem2
  # library(ggnewscale)

  id1 = idtab %>% filter(group ==unique(data3$group)[i]) %>%.[1:(floor(nrow(data)/2)),] %>%.$label
  id2 = idtab %>% filter(group ==unique(data3$group)[i]) %>%.[(floor(nrow(data)/2) + 1):nrow(data),] %>%.$label


  p = p +
    geom_curve(data = tab2 %>% filter(label %in% id1),
               aes_string(x = "xend",
                          y = "yend",
                          xend = "x",
                          yend = "y",
                          group = "groupbc",
                          color ="groupr",
                          size = "size"),
              curvature = -corva
               ) +
    geom_curve(data = tab2 %>% filter(label %in% id2),
               aes_string(x = "xend",
                          y = "yend",
                          xend = "x",
                          yend = "y",
                          group = "groupbc",
                          color ="groupr",
                          size = "size"),
               curvature = corva
    ) +
    scale_color_manual(values = c("#91331FFF","#46732EFF")) +
    # guides(shape = guide_legend(override.aes = list(fill = "blue") )) +
    ggnewscale::new_scale_fill() +
    geom_point(data = tab2,
                   aes_string(x = "x",
                              y = "y",
                  color = "groupr",
                  fill = "groupr",
                  size = tem2*0.65
                  ),
                  pch = 21) +
    scale_color_manual(values = c("#91331FFF","#46732EFF")) +
    scale_fill_manual(values = c("#91331FFF","#46732EFF")) +
    ggnewscale::new_scale_fill() +
    geom_point(data = topdat,aes(x = xend, y = yend),pch = numpoint2,size =8,
               color = "black",fill = "#FFF5EB")

  p

}

# env.dat = env
# tabOTU = tabOTU1
# MetalTast(env.dat = env,tabOTU = tabOTU1,distance = "jcd",method = "metal")

