






cor_link.dat <- function(data,
                      envdata = report,
                      Ptab,
                      sig = F,
                      zoom = 2,
                      show = "z",
                      corva = 0,
                      range = 1,
                      lacx = "left",
                      lacy = "bottom",
                      numpoint2 = 21,
                      curvature = 0.2,
                      p.thur = 0.3,
                      onlysig = TRUE
){


  if (show == "y"&lacx == "right") {
    data = data.frame(x = max(data$x)+1,y = data$y,label = data$label)

  }else if(show == "y"&lacx == "left"& lacy == "bottom"){
    data = data.frame(x = min(data$x)-1,y = data$y-1,label = data$label)
  }else if(show == "y"&lacx == "left"& lacy == "top"){
    data = data.frame(x = min(data$x)-1,y = data$y,label = data$label)
  }else if(show == "x"& lacy == "top"&lacx == "right"){
    data = data.frame(x = data$x ,y =max(data$y)+1,label = data$label)
  }else if(show == "x"& lacy == "top"&lacx == "left"){
    data = data.frame(x = data$x ,y =max(data$y)+1,label = data$label)
  }else if(show == "x"& lacy == "bottom"&lacx == "left"){
    data = data.frame(x = data$x,y =min(data$x)-1,label = data$label)
  }else if(show == "x"& lacy == "bottom"&lacx == "right"){
    data = data.frame(x = data$x ,y =min(data$x)-1,label = data$label)
  } else if(show == "z"){
    data = data
  }



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
  if (onlysig == TRUE& sig) {
    data3 = data3 %>% filter(p < p.thur )
  }

  # ggplot() + geom_point(aes(x,y),data = data) +
  #   geom_point(aes(x -1,y +1),data = data)
  #
  # ggplot() + geom_point(aes(x,y),data = data) +
  #   geom_point(aes(xend,yend),data = topdat)
  n = x +1

  if (show == "y"&lacx == "left") {
    dat = data.frame(x =  min(data$x) - zoom,y = data$y )
    seqnum = (max(dat$y) - min(dat$y))/n
    tem = seq(from=min(dat$y), to=max(dat$y),by=seqnum)
    tem = tem[length(tem):1]

    topdat = data.frame(x=seq(from=min(dat$x), to=max(dat$x),by=seqnum),
                        y = tem) %>% slice(-1)
    topdat = topdat[-nrow(topdat),]
    colnames(topdat) = paste0(colnames(topdat),"end")
    topdat$group = unique(data3$group)

  } else if (show == "y"&lacx == "right"){
    dat = data.frame(x = max(data$x) + zoom,y = data$y )
    seqnum = (max(dat$y) - min(dat$y))/n
    topdat = data.frame(x= max(data$x) + zoom,
                        y = seq(from=min(dat$y), to=max(dat$y),by=seqnum)) %>%
      dplyr::slice(-1)
    topdat = topdat[-nrow(topdat),]
    colnames(topdat) = paste0(colnames(topdat),"end")
    topdat$group = unique(data3$group)
  } else if(show == "x"&lacy == "bottom"){
    dat = data.frame(x = data$x,y = min(data$y) - zoom)
    seqnum = (max(dat$x) - min(dat$x))/n
    tem = seq(from=min(dat$y), to=max(dat$y),by=seqnum)
    tem = tem[length(tem):1]

    topdat = data.frame(x=seq(from=min(dat$x), to=max(dat$x),by=seqnum),
                        y = tem) %>% slice(-1)
    topdat = topdat[-nrow(topdat),]
    colnames(topdat) = paste0(colnames(topdat),"end")
    topdat$group = unique(data3$group)

  } else if( show == "x"&lacy == "top"){
    dat = data.frame(x = data$x ,y = max(data$y)+zoom)
    seqnum = (max(dat$x) - min(dat$x))/n
    tem = seq(from=min(dat$y), to=max(dat$y),by=seqnum)
    # tem = tem[length(tem):1]

    topdat = data.frame(x=seq(from=min(dat$x), to=max(dat$x),by=seqnum),
                        y = tem) %>% slice(-1)
    topdat = topdat[-nrow(topdat),]
    colnames(topdat) = paste0(colnames(topdat),"end")
    topdat$group = unique(data3$group)

  }
  if (show=="z"&lacx == "left"& lacy == "bottom") {
    dat = data.frame(x = data$x + zoom,y = data$y + zoom)
    seqnum = (max(dat$x) - min(dat$x))/n
    tem = seq(from=min(dat$y), to=max(dat$y),by=seqnum)
    tem = tem[length(tem):1]

    topdat = data.frame(x=seq(from=min(dat$x), to=max(dat$x),by=seqnum),
                        y = tem) %>% slice(-1)
    topdat = topdat[-nrow(topdat),]
    colnames(topdat) = paste0(colnames(topdat),"end")
    topdat$group = unique(data3$group)

  } else if (show=="z"&lacx == "right"& lacy == "bottom"){
    dat = data.frame(x = data$x - zoom,y = data$y + zoom)
    seqnum = (max(dat$x) - min(dat$x))/n
    topdat = data.frame(x=seq(from=min(dat$x), to=max(dat$x),by=seqnum),
                        y = seq(from=min(dat$y), to=max(dat$y),by=seqnum)) %>%
      dplyr::slice(-1)
    topdat = topdat[-nrow(topdat),]
    colnames(topdat) = paste0(colnames(topdat),"end")
    topdat$group = unique(data3$group)
  } else if(show=="z"&lacx == "right"&lacy == "top"){
    dat = data.frame(x = data$x - zoom,y = data$y - zoom)
    seqnum = (max(dat$x) - min(dat$x))/n
    tem = seq(from=min(dat$y), to=max(dat$y),by=seqnum)
    tem = tem[length(tem):1]

    topdat = data.frame(x=seq(from=min(dat$x), to=max(dat$x),by=seqnum),
                        y = tem) %>% slice(-1)
    topdat = topdat[-nrow(topdat),]
    colnames(topdat) = paste0(colnames(topdat),"end")
    topdat$group = unique(data3$group)

  } else if( show=="z"&lacx == "left"&lacy == "top"){
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

  id1 = idtab %>% filter(group ==unique(data3$group)[1]) %>%.[1:(floor(nrow(data)/2)),] %>%.$label
  id2 = idtab %>% filter(group ==unique(data3$group)[1]) %>%.[(floor(nrow(data)/2) + 1):nrow(data),] %>%.$label
  print(id1)
  return(list(tab2,id1,id2))

}



