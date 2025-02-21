
cor_linkj = function(
    data,
    p,
    repR,
    repP,
    matrix.line = list(),
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

  if (lacx == "right"&lacy == "bottom") {
    corva = -corva
  }

  id.x = matrix.line$x

  if (!is.null(id.x)) {
    repR.x = repR[,c("Envs",paste0(id.x,"r.BC"))]
    repP.x = repP[,c("Envs",paste0(id.x,"p.BC"))]

    dat = cor_link.dat (data = data,
                        envdata = repR.x,
                        Ptab= repP.x,
                        show = "x",
                        zoom = zoom,
                        corva = corva,
                        range = range,
                        lacx = lacx,
                        lacy = lacy,
                        numpoint2 = 21,
                        curvature = 0.2,
                        p.thur = 0.3,
                        onlysig = TRUE
    )
    dat = dat[[1]]
    dat$class = "x"
    tem2.0 = dat$size
    topdat = data.frame(xend = unique(dat$xend),yend = unique(dat$yend),group = unique(dat$group))


  }else{
    dat = NULL
    tem2.0 = NULL
    topdat = NULL
  }


  id.y = matrix.line$y
  if (!is.null(id.y)) {
    # repR.y = repR[,c("Envs",colnames(repR)[str_detect(colnames(repR),id.y)])]
    # repP.y = repP[,c("Envs",colnames(repP)[str_detect(colnames(repP),id.y)])]
    repR.y = repR[,c("Envs",paste0(id.y,"r.BC"))]
    repP.y = repP[,c("Envs",paste0(id.y,"p.BC"))]


    dat1 = cor_link.dat (data = data,
                         envdata = repR.y,
                         Ptab= repP.y,
                         show = "y",
                         zoom = zoom,
                         corva = -corva,
                         range = range,
                         lacx = lacx,
                         lacy = lacy,
                         numpoint2 = 21,
                         curvature = 0.2,
                         p.thur = 0.3,
                         onlysig = TRUE
    )
    dat1 = dat1[[1]]
    dat1$class = "y"
    tem2.1 = dat1$size
    topdat1 = data.frame(xend = unique(dat1$xend),yend = unique(dat1$yend),group = unique(dat1$group))

  } else {
    dat1 = NULL
    tem2.1 = NULL
    topdat1 = NULL
  }



  id.z = matrix.line$z
  if (!is.null(id.z)) {
    # repR.z = repR[,c("Envs",colnames(repR)[str_detect(colnames(repR),id.z)])]
    # repP.z = repP[,c("Envs",colnames(repP)[str_detect(colnames(repP),id.z)])]
    repR.z = repR[,c("Envs",paste0(id.z,"r.BC"))]
    repP.z = repP[,c("Envs",paste0(id.z,"p.BC"))]

    dat20 = cor_link.dat (data = data,
                          envdata = repR.z,
                          Ptab= repP.z,
                          show = "z",
                          zoom = zoom,
                          corva = -0.2,
                          range = range,
                          lacx = lacx,
                          lacy = lacy,
                          numpoint2 = 21,
                          curvature = 0.2,
                          p.thur = 0.3,
                          onlysig = TRUE
    )
    dat2 = dat20[[1]]
    dat2$class = "z"
    tem2.2 = dat2$size
    topdat2 = data.frame(xend = unique(dat2$xend),yend = unique(dat2$yend),group = unique(dat2$group))


  }else{
    dat2 = NULL
    tem2.2 = NULL
    topdat2 = NULL
  }



  dat0 = rbind(dat,dat1,dat2)

  tem2 = c(tem2.0,tem2.1,tem2.2)
  head(dat0)

  topdat0 = rbind(topdat ,topdat1 ,topdat2 )

  datxy = dat0 %>% dplyr::filter(class %in%c("x","y"))

  id1 = dat20[[2]]
  id2 = dat20[[3]]

  if (is.null(dat2)) {

  }



  if (!is.null(datxy)) {
    p1 = p+
      # new_scale("size") +
      geom_curve(data = datxy  %>% filter(label %in% id1),
                 aes_string(x = "xend",
                            y = "yend",
                            xend = "x",
                            yend = "y",
                            group = "groupbc",
                            color ="groupr",
                            size = "size"),
                 curvature = corva,show.legend = TRUE
      )

  }else{
    p1 = p
  }
  if (!is.null(dat2)) {
    p2 = p1 +
      geom_curve(data = dat2  %>% filter(label %in% id1),
                 aes_string(x = "xend",
                            y = "yend",
                            xend = "x",
                            yend = "y",
                            group = "groupbc",
                            color ="groupr",
                            size = "size"),
                 curvature = -corva,show.legend = TRUE
      )
  }else{
    p2 = p1
  }


  if (!is.null(datxy)) {
    p3 = p2+ ggnewscale::new_scale("size") +

      geom_curve(data = datxy  %>% filter(label %in% id2),
                 aes_string(x = "xend",
                            y = "yend",
                            xend = "x",
                            yend = "y",
                            group = "groupbc",
                            color ="groupr",
                            size = "size"),
                 curvature = -corva,show.legend = FALSE
      )
  }else{
    p3 = p2
  }
  if (!is.null(dat2)) {
    p4 = p3 + geom_curve(data = dat2  %>% filter(label %in% id2),
                         aes_string(x = "xend",
                                    y = "yend",
                                    xend = "x",
                                    yend = "y",
                                    group = "groupbc",
                                    color ="groupr",
                                    size = "size"),
                         curvature = corva,show.legend = FALSE
    )

  } else{
    p4 = p3
  }


  p5 = p4 +
    scale_size_continuous(name = "Size.matel") +
    scale_color_manual(values = c("#91331FFF","#46732EFF")) +
    # guides(shape = guide_legend(override.aes = list(fill = "blue") )) +
    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +
    geom_point(data = dat0,
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
    geom_point(data = topdat0,aes(x = xend, y = yend),pch = numpoint2,size =8,
               color = "grey60",fill = "#FFF5EB") +
    ggrepel::geom_text_repel(data = topdat0,aes(x = xend, y = yend,label = group))


  return(p5)
}








