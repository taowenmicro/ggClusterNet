#' Combine microbiome data by classification level
#'
#' @param ps phyloseq object
#' @param  ranks tax of microbiome data,could be one of Phylum, Order, Class, Genus ta al.
#' @examples
#' tax_glom_wtTop(ps = ps,ranks = "Genus",Top = 10)
#' tax_glom_wtTop(ps = ps,ranks = "Class",Top = 10)
#' tax_glom_wtTop(ps = ps,ranks = "Order",Top = 10)
#' tax_glom_wtTop(ps = ps,ranks = "Phylum",Top = 10)
#' @return S4 phyloseq abject
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export


tax_glom_wtTop = function(ps = ps,
                       ranks = "Phylum",
                       Top = 10
                       ){
    if (  is.numeric(ranks)) {
      ranks <- phyloseq::rank.names(ps)[ranks]
    }


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

    if (is.vector(taxcon)) {
      taxcon = data.frame(row.names = taxcon,ranks = taxcon)
      colnames(taxcon) = ranks
    }

    #-tax name with NA wound be repeated with unknown
    taxcon[[ranks]][is.na(taxcon[[ranks]])] = "unknown"
    row.names(taxcon) <- taxcon[[ranks]]
    head(otucon)
    # tem = rowSums(otucon) %>%
    #   as.data.frame()
    # colnames(tem) = "count"
    # id = tem %>% arrange(desc(count)) %>% head(Top) %>% row.names()
    # head(taxcon)


    pscon <- phyloseq::phyloseq(
      phyloseq::otu_table( as.matrix(otucon),taxa_are_rows = TRUE),
      phyloseq::tax_table(as.matrix(taxcon)),
      phyloseq::sample_data(ps)
    )

    otu = phyloseq::otu_table(pscon)
    tax = phyloseq::tax_table(pscon)
    otu[is.na(otu)] = 0
    i = 2
    for (i in 1:dim(tax)[1]) {
      if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
        tax[i,ranks] =tax[i,ranks]
      } else {
        tax[i,ranks]= "Others"
      }
    }
    phyloseq::tax_table(pscon)= tax

    pscon2 = pscon %>%tax_glom_wt(ranks = ranks)

    # pscon2 %>% vegan_tax()

    return(pscon2)
  }

