
matcorplotj = function(

  env.dat = env1,
  tabOTU = tabOTU1,
  corva = 0.05,
  zoom = 3,
  diag = F,
  range = 0.01,
  numpoint = 21,
  numpoint2 = 20,
  sig = F,
  siglabel = FALSE,
  shownum = F,
  curvature = 0,
  numsymbol = NULL,
  lacx = "left",
  lacy = "top",
  p.thur = 0.05,
  onlysig = TRUE,
  method.cor = "spearman",
  method = "metal",
  distance = "bray",
  cor.p = 0.05,
  x = TRUE,
  y = TRUE,
  matrix.line = list(x = c(names(tabOTU1)[1:2]),
                     y = c(names(tabOTU1)[3]),
                     z = c(names(tabOTU1)[3]))

){

  rep = MetalTast(env.dat = env.dat, tabOTU = tabOTU,distance = distance,method = method)
  repR = rep[c(-seq(from=1,to=dim(rep)[2],by=2)[-1])]
  repP = rep[seq(from=1,to=dim(rep)[2],by=2)]

  result <- Miccorplot3(data = env.dat,
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

  cor_linkj(data = result[[2]],# env table
            p = result[[1]], # ggplot object
            repR = repR,# R value table
            repP = repP,#p value table
            zoom = zoom,
            matrix.line = matrix.line,
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

