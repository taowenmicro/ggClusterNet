
#' Combine microbiome data by classification level
#'
#' @param ps phyloseq object
#' @param  ranks tax of microbiome data,could be one of Phylum, Order, Class, Genus ta al.
#' @examples
#' tax_glom_wt(ps = ps,ranks = "Genus")
#' tax_glom_wt(ps = ps,ranks = "Class")
#' tax_glom_wt(ps = ps,ranks = "Order")
#' tax_glom_wt(ps = ps,ranks = "Phylum")
#' @return S4 phyloseq abject
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


tax_glom_wt <- function(ps = ps,ranks = "Phylum") {


  otu <- as.data.frame(t(vegan_otu(ps)))
  tax <- as.data.frame(vegan_tax(ps))

  # building group
  tax[[ranks]][is.na(tax[[ranks]])] = "Unknown"
  tax[[ranks]][tax[[ranks]] == ""] = "Unknown"
  tax[[ranks]][tax[[ranks]] == "NA"] = "Unknown"
  split <- split(otu,tax[[ranks]])
  #calculate sum by group
  apply <- lapply(split,function(x)colSums(x[]))
  # result chack
  otucon <- do.call(rbind,apply)

  taxcon <- tax[1:match(ranks,colnames(tax))]
  taxcon <- taxcon[!duplicated(tax[[ranks]]),]
  #-tax name with NA wound be repeated with unknown
  taxcon[[ranks]][is.na(taxcon[[ranks]])] = "unknown"
  row.names(taxcon) <- taxcon[[ranks]]


  pscon <- phyloseq(
    otu_table( as.matrix(otucon),taxa_are_rows = TRUE),
    tax_table(as.matrix(taxcon)),
    sample_data(ps)
  )
  return(pscon)
}
