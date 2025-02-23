
# otupath = "../ggClusterNet_Decument/2025.02.19-ggClusterNet2文章梳理/"
# fs::dir_create(otupath)
# netpath = paste(otupath,"/result.micro.cor.other/",sep = "")
# dir.create(netpath)
#
# otu1 = as.data.frame(t(ggClusterNet::vegan_otu(ps)))
# head(otu1)
# # 无论怎么绘制，都这么计算关系
# tabOTU1 = list(bac = otu1,
#                fun = otu1,
#                pre = otu1
# )
#
#
#
#
# rep = MetalTast(env.dat = env1, tabOTU = tabOTU1,
#                 distance = "bray",
#                 method = "metal"
#
#
# )
# repR = rep[c(-seq(from=1,to=dim(rep)[2],by=2)[-1])]
# repP = rep[seq(from=1,to=dim(rep)[2],by=2)]
#
#
# #   验证方向的函数#---------
# show = "y"
#
# res = two.cor (
#   env.dat1 = env.dat,
#   env.dat2 = env.dat,
#   lacx1 = "left",
#   lacy1 = "top",
#   lacx2 = "right",
#   lacy2 = "bottom",
#   show = show,
#   zoom.x = 15,
#   zoom.y = 15
# )
#
#
# # data = res[[2]]# env table
# # p = res[[1]] # ggplot object
# # envdata = repR# R value table
# # Ptab = repP#p value table
# # zoom = 3
# # sig = F
# # show = show
# # topdat = NULL
# # corva = 0
# # numpoint2 = 21
# # curvature = 0
# # range = 0.01# line size
# # res = res
# # p.thur = 0.05#sig value
# # onlysig = F # if show sig connect
#
# cor_link.two(
#   data = res[[2]],# env table
#   p = res[[1]], # ggplot object
#   envdata = repR,# R value table
#   Ptab = repP,#p value table
#   zoom = 3,
#   sig = F,
#   show = show,
#   topdat = NULL,
#   corva = 0,
#   numpoint2 = 21,
#   curvature = 0,
#   range = 0.01,# line size
#   res = res,
#   p.thur = 0.05,#sig value
#   onlysig = F # if show sig connect
#
#
# )
cor_link.two = function(
    data = res[[2]],# env table
    p = res[[1]], # ggplot object
    envdata = repR,# R value table
    Ptab = repP,#p value table
    zoom = 3,
    sig = F,
    show = "x",
    topdat = NULL,
    corva = 0,
    numpoint2 = 21,
    curvature = 0,
    range = 0.01,# line size
    res = res,
    p.thur = 0.05,#sig value
    onlysig = F # if show sig connect



){

  p = cor_linkx(
    data = res[[2]],# env table
    p = res[[1]], # ggplot object
    envdata = repR,# R value table
    Ptab = repP,#p value table
    zoom = zoom,
    show = show,
    topdat = NULL,
    corva = corva,
    numpoint2 = numpoint2,
    curvature = curvature,
    range = range,# line size
    lacx = res[[4]],#  left or right
    lacy = res[[6]],# top or bottom
    p.thur = p.thur,#sig value
    onlysig = onlysig # if show sig connect
  )

  p[[1]]

  p1 = cor_linkx(
    data = res[[3]],# env table
    p = res[[1]], # ggplot object
    envdata = repR,# R value table
    Ptab = repP,#p value table
    zoom = zoom,
    show = show,
    topdat = NULL,
    corva = corva,
    numpoint2 = numpoint2,
    curvature = curvature,
    range = range,# line size
    lacx = res[[5]],#  left or right
    lacy = res[[7]],# top or bottom
    p.thur = p.thur,#sig value
    onlysig = onlysig # if show sig connect
  )
  p1[[1]]

  topdat = p[[2]]
  topdat1 = p1[[2]]


  if (show == "y") {
    a = c(topdat$xend +topdat1$xend[length(topdat1$xend):1])/2
    b = c(topdat$yend +topdat1$yend[length(topdat1$yend):1])/2
    topdat$xend = a
    topdat$yend = b
  } else{
    a = c(topdat$xend +topdat1$xend)/2
    b = c(topdat$yend +topdat1$yend)/2
    topdat$xend = a
    topdat$yend = b
  }




  p2 = cor_linkx(
    data = res[[3]],# env table
    p = res[[1]], # ggplot object
    envdata = repR,# R value table
    Ptab = repP,#p value table
    topdat = topdat,
    zoom = zoom,
    show = show,
    corva = corva,
    numpoint2 = numpoint2 ,
    curvature = curvature,
    range = range,# line size
    lacx = res[[5]],#  left or right
    lacy = res[[7]],# top or bottom
    p.thur = p.thur,#sig value
    onlysig = onlysig # if show sig connect
  )
  p2[[1]]


  p3 = cor_linkx(
    data = res[[2]],# env table
    p = p2[[1]], # ggplot object
    envdata = repR,# R value table
    Ptab = repP,#p value table
    topdat = p2[[2]],
    zoom = zoom,
    show = show,
    corva = corva,
    numpoint2 = numpoint2,
    curvature = curvature,
    range = range,# line size
    lacx = res[[4]],#  left or right
    lacy = res[[6]],# top or bottom
    p.thur = p.thur,#sig value
    onlysig = onlysig # if show sig connect
  )


  # p3[[1]]
  return(p3[[1]])
}
