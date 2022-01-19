#' Stacked histograms add error bars, marked by significance
#'
#' @param ps a object of phyloseq
#' @param  j rank of tax
#' @param Top top number of number selected
#' @param errbar add a distinctive label, TRUE or FEASE could be selected
#' @param result output from aovMcomper or KwWlx. You can also import result calculated from other software (a data frame)
#' @param add_abc add sig abc T or F
#' @examples
#' ps_merge_Tax (ps = ps,j = "Phylum",Top = 10)
#' @return list
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn} PengHao Xie \email{2019103106@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export
#

ps_merge_Tax = function(ps = ps,
                        j = "Phylum",
                        Top = 10
){
  ps.tem = ps %>% tax_glom_wt(ranks = j)
  otu = ps.tem %>% vegan_otu() %>% t() %>% as.data.frame()
  tax = ps.tem %>% vegan_tax() %>% as.data.frame()
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {

      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  head(tax)
  tax$Phylum %>% unique()
  tax_table(ps.tem)= as.matrix(tax)
  ps.j = ps.tem %>% tax_glom_wt(ranks = j)
  return(ps.j)
}
