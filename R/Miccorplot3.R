

Miccorplot3 <- function(data,
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
    geom_tile(aes(x = x, y = y),data = df0,fill=NA,color='gray',linewidth =0.5) +
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
    if (lacy == "bottom"&lacx == "right") {
      labdat$x = length(labdat$x) - labdat$x
      labdat$z = length(labdat$x) - labdat$z
      p <- p +  geom_text(aes(x = x+1 , y = z +1,label=labx),labdat,
                          vjust =  1,hjust = 0.5)+
        theme(plot.margin=unit(rep(2,4),'cm'))
      lindat = data.frame(x = labdat$x + 1,y = labdat$z +1 ,label = labdat$labx)
    } else if(lacy == "top"&lacx == "right"){
      # labdat$x = length(labdat$x) - labdat$x
      # labdat$z = length(labdat$x) - labdat$z
      p <- p +  geom_text(aes(x = x+1 , y = z -1,label=labx),labdat,
                          vjust =  1,hjust = 0.5)+
        theme(plot.margin=unit(rep(2,4),'cm'))
      lindat = data.frame(x = labdat$x + 1,y = labdat$z - 1 ,label = labdat$labx)
    }  else if(lacy == "bottom"&lacx == "left") {
      labdat$x = length(labdat$x) - labdat$x
      labdat$z = length(labdat$x) - labdat$z
      p <- p +  geom_text(aes(x = x +1, y = z+1 ,label=labx),labdat,
                          vjust =  1,hjust = 0.5)+
        theme(plot.margin=unit(rep(2,4),'cm'))
      lindat = data.frame(x = labdat$x + 1,y = labdat$z +1 ,label = labdat$labx)
    } else if(lacy == "top"&lacx == "left"){
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

