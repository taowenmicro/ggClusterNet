
two.cor2 = function(
    env.dat1 = env.dat,
    env.dat2 = env.dat,
    numpoint = numpoint,
    numsymbol = NULL,
    lacx1 = "left",
    lacy1 = "bottom",
    lacx2 = "right",
    lacy2 = "top",
    show = "x",
    zoom.x = 8,
    zoom.y = 8
){

  if (show == "x") {
    lacy1 = "top"
    zoom.x = 0
    zoom.y = abs(zoom.y)
    lacy2 = "bottom"
    # if (lacx2 == "right") {
    #
    # }
  } else if(show == "y"){
    lacx1 = "right"
    zoom.y = 0
    zoom.x = abs(zoom.x)
    lacx2 = "left"
  } else if (show == "z"){
    lacy1 = "bottom"
    lacy2 = "top"
    # zoom.x = abs(zoom.x)
    # zoom.y = abs(zoom.y)
  }

  result <- Miccorplot3(data = env.dat1,
                        method.cor = "spearman",
                        cor.p = 0.05,
                        x = TRUE,
                        y = TRUE,
                        diag = TRUE,
                        lacx = lacx1,
                        lacy = lacy1,
                        sig = TRUE,
                        siglabel = FALSE,
                        shownum = F,
                        numpoint = numpoint,
                        numsymbol = numsymbol)

  result[[1]]
  link1 = result[[2]]
  #  当前这两个表格需要平移
  da1 = result[[3]]
  labdat = result[[4]]


  result <- Miccorplot3(data = env.dat2,
                        method.cor = "spearman",
                        cor.p = 0.05,
                        x = TRUE,
                        y = TRUE,
                        diag = TRUE,
                        lacx = lacx2,
                        lacy = lacy2,
                        sig = TRUE,
                        siglabel = FALSE,
                        shownum = F,
                        numpoint = numpoint,
                        numsymbol = numsymbol)

  result[[1]]
  link2 = result[[2]]
  #  当前这两个表格需要平移
  da2 = result[[3]]
  labdat = result[[4]]
  if (show == "z"&lacx1 == "left"&lacx2 == "right") {
    zoom.y = abs(zoom.y)
    zoom.x = abs(zoom.x)
    da2$x = da2$x + zoom.x
    da2$y = da2$y + zoom.y

  } else if (show == "z"&lacx1 == "left"&lacx2 == "left") {
    zoom.y = abs(zoom.y)
    zoom.x = abs(zoom.x)
    da2$x = da2$x
    da2$y = da2$y + zoom.y

  }else if (show == "z"&lacx1 == "right"&lacx2 == "left") {
    zoom.y = abs(zoom.y)
    zoom.x = -abs(zoom.x)

    da2$x = da2$x + zoom.x
    da2$y = da2$y + zoom.y
  }else if (show == "z"&lacx1 == "right"&lacx2 == "right") {
    zoom.y = abs(zoom.y)
    zoom.x = 0

    da2$x = da2$x
    da2$y = da2$y + zoom.y
  }

  if (show == "x"&lacx2 == "left") {
    da2$x = da2$x
    da2$y = da2$y + zoom.y
  }else if (show == "x"&lacx2 == "right"){
    da2$x = da2$x
    da2$y = da2$y + zoom.y
  }else if (show == "y"){
    da2$x = da2$x + zoom.x
    da2$y = da2$y
  }else if(show == "x"&lacy2 == "right"){
    da2$x = da2$x
    da2$y = da2$y + zoom.y
  }

  p1 = ggplot(data = rbind(da1,da2)) +
    geom_tile(aes(x = x, y = y),fill=NA,color='gray',linewidth =0.5) +
    scale_size(range = c(0.1, 6))+
    scale_fill_distiller(palette="RdYlBu") + theme_void() +
    geom_point(aes(x = x, y = y,
                   size=r ,fill=r), shape=numpoint,color='white')+
    scale_fill_gradient2(low = "#008B45FF",
                         high = "#EE0000FF") +
    theme(plot.margin=unit(rep(2,4),'cm'))
  p1
  link2$x = link2$x + zoom.x
  link2$y = link2$y + zoom.y

  return(list(p1,link1,link2,lacx1,lacx2,lacy1,lacy2))
}



