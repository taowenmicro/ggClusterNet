#' Construct edge files, add weights, positive and negative correlation and other basic information
#'
#' @param cor Correlation matrix
#' @param node Node file, containing calculated node coordinates
#' @examples
#' data(ps)
#' result = corMicro (ps = ps,N = 120,r.threshold=0.8,p.threshold=0.05,method = "pearson")
#' #Extract correlation matrix
#' cor = result[[1]]
#' # building the node group
#' netClu = data.frame(ID = row.names(cor),group =rep(1:3,length(row.names(cor)))[1:length(row.names(cor))] )
#' netClu$group = as.factor(netClu$group)
#' #Calculate node location
#' result2 = PolygonRrClusterG(cor = cor,nodeGroup = netClu )
#' node = result2[[1]]
#' edge = edgeBuild(cor = cor,plotcord = node)
#'
#' @return edge which contains OTU and its coordinates, the correlation of edge between nodes
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} Penghao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


edgeBuild = function(cor = cor,node = node){
  tem1 = cor %>%
    tidyfst::mat_df() %>%
    dplyr::filter(row != col) %>%

    dplyr::rename(OTU_1 = row,OTU_2 = col,weight = value ) %>%
    dplyr::filter(weight != 0)
  head(tem1)

  tem2 = tem1 %>% dplyr::left_join(node,by = c("OTU_1" = "elements")) %>%
    dplyr::rename(Y1 = X2)
  head(tem2)
  tem3 = node %>%
    dplyr::rename(Y2 = X2,X2 = X1) %>%
    dplyr::right_join(tem2,by = c("elements" = "OTU_2")) %>%
    dplyr::rename(OTU_2 = elements)

  edge = tem3 %>%
    dplyr::mutate(
      cor = ifelse(weight > 0,"+","-")
    )
  colnames(edge)[8] = "cor"

  return(edge)
}
