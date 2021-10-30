#' Use microbiome data to map petals
#'
#' @param otu otu table of microbiome. data.frame
#' @param tax taxonmy table of microbiome. data.frame
#' @param map table.data.frame
#' @param ps phyloseq Object, contains OTU tables, tax table and map table, represented sequences,phylogenetic tree.
#' @param group colnames which selected for show
#' @param rep repeat number of microbial data.
#' @param m1 Petal shape, square to round to prismatic, the value gradually decreases
#' @param start The rotation angle of the petals, the greater the value, the greater the angle
#' @param a The width of the petals
#' @param b Distance from petal to center
#' @param lab.leaf The distance from the label to the center of the circle
#' @param col.cir Center color
#' @examples
#' library(phyloseq)
#' library(ggplot2)
#' data(ps)
#'
#' map = as.data.frame(sample_data(ps))
#' map$Group1 <- c("A","B","C","D","E","F ")
#' sample_data(ps) = map
#'
#' p <- ggflower(ps = ps,
#'               rep = 3,
#'               group = "Group1",
#'               start = 1,
#'               m1 = 2,
#'               a = 0.3,
#'               b = 0.3,
#'               lab.leaf = 1,
#'               col.cir = "yellow",
#'               b.cir = 0.8,
#'               a.cir = 0.8
#' )
#' p2 <- p + scale_fill_brewer(palette = "Paired")
#' p2
#' p <- ggflower(ps = ps,
#'               rep = 3,
#'               group = "Group1",
#'               start = 1,
#'               m1 = 1,
#'               a = 0.3,
#'               b = 1,
#'               lab.leaf = 1,
#'               col.cir = "yellow"
#' )
#' p3 <- p + scale_fill_brewer(palette = "Paired")
#'
#' p3
#' p5 <- ggflower(ps = ps,
#'                rep = 3,
#'                group = "Group1",
#'                start = 1,
#'                m1 = 2,
#'                a = 0.2,
#'                b = 1,
#'                lab.leaf = 1,
#'                col.cir = "yellow"
#' )
#' p5
#'
#' p <- ggflower(ps = ps,
#'               rep = 3,
#'               group = "Group1",
#'               start = 30,
#'               m1 = 2,
#'               a = 0.2,
#'               b = 1,
#'               lab.leaf = 1,
#'               col.cir = "yellow"
#' )
#' p1 <- p + scale_fill_brewer(palette = "Paired")
#'
#' p1
#' @return ggplot objects
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


ggflower = function(otu = NULL,tax = NULL,map = NULL,ps = NULL,
                    group = "Group",
                    rep = 6,
                    m1 = 2,
                    start = 1,
                    a = 0.2,
                    b = 1,
                    lab.leaf = 1,
                    col.cir = "yellow",
                    a.cir = 0.5,
                    b.cir = 0.5,
                    m1.cir = 2
                    ){
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  mapping = as.data.frame(sample_data(ps))

  aa = vegan_otu(ps)
  otu_table = as.data.frame(t(aa))
  count = aa
  sub_design <- as.data.frame(sample_data(ps))


  sub_design$SampleType = sub_design$Group
  sample_data(ps ) = sub_design


  #
  pick_val_num <- rep *2/3
  count[count > 0] <- 1
  count2 = as.data.frame(count)


  #dive group
  iris.split <- split(count2,as.factor(sub_design$Group))
  # cul mean
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  # conbine result
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine)

  ven2[ven2 < pick_val_num]  = 0
  ven2[ven2 >=pick_val_num]  = 1
  ven2 = as.data.frame(ven2)


  #
  ven3 = as.list(ven2)
  ven2 = as.data.frame(ven2)

  all_num = dim(ven2[rowSums(ven2) == length(levels(sub_design$Group)),])[1]

  ven2[,1] == 1

  A = rep("A",length(colnames(ven2)))
  B = rep(1,length(colnames(ven2)))

  i = 1
  for (i in 1:length(colnames(ven2))) {
    B[i] = length(ven2[rowSums(ven2) == 1,][,i][ven2[rowSums(ven2) == 1,][,i] == 1])
    A[i] = colnames(ven2)[i]
  }

  n   <- length(A)
  deg <- 360 / n
  t = 1:n
  p <- ggplot() +
    # geom_point(aes(x = 5 + cos((start + deg * (t - 1)) * pi / 180) * lab.leaf, y = 5 + sin((start + deg * (t - 1)) * pi / 180) *lab.leaf)) +
    ggforce::geom_ellipse(aes(x0 = 5 + cos((start + deg * (t - 1)) * pi / 180),
                     y0 = 5 + sin((start + deg * (t - 1)) * pi / 180),
                     a = a,
                     b = b,
                     angle = (n/2 +seq(0,1000,2)[1:n])/n * pi,
                     m1 = m1,
                     fill = as.factor(1:n)),show.legend = F) +
    ggforce::geom_ellipse(aes(x0 = 5,y0 = 5,a = a.cir,b = b.cir,angle = 0,m1 = m1.cir),fill = col.cir) +
    geom_text(aes(x = 5,y = 5,label = paste("OVER :",all_num,sep = ""))) +
    geom_text(aes(
      x = 5 + cos((start + deg * (t - 1)) * pi / 180) * lab.leaf,
      y = 5 + sin((start + deg * (t - 1)) * pi / 180) * lab.leaf,
      label = paste(A,":",B,sep = "")),angle = 360/n*((1:n)-1)  ) +
    coord_fixed() + theme_void()
p
return(p)
}




