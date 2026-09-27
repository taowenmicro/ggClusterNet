
#' Fast module roles
#'
#' Uses an existing vertex attribute "module" when membership is not supplied.
#' @export
module.roles <- function(
    comm_graph,
    membership = NULL,
    role_weights = c("degree", "strength"),
    zi_cut = 2.5,
    pi_cut = 0.62
) {
  role_weights <- match.arg(role_weights)

  if (is.null(membership)) {
    membership <- igraph::vertex_attr(comm_graph, "module")
    if (is.null(membership))
      stop("No module membership supplied and vertex attribute 'module' is absent.")
  }

  .zipi_roles_fast(
    comm_graph,
    membership = membership,
    role_weights = role_weights,
    zi_cut = zi_cut,
    pi_cut = pi_cut
  )
}



.zipi_roles_fast <- function(
    g,
    membership,
    role_weights = c("degree", "strength"),
    zi_cut = 2.5,
    pi_cut = 0.62
) {
  role_weights <- match.arg(role_weights)

  if (igraph::is_directed(g))
    g <- igraph::as_undirected(g, mode = "collapse")

  n <- igraph::vcount(g)
  m <- igraph::ecount(g)
  taxa <- igraph::V(g)$name
  if (is.null(taxa)) taxa <- as.character(seq_len(n))

  module <- match(as.character(membership), unique(as.character(membership)))
  if (length(module) != n) stop("membership length mismatch.")

  ## No-edge graph.
  if (m == 0L) {
    out <- data.frame(
      z = 0, module = module, p = 0,
      roles = "Peripherals", role_7 = "ultra peripheral",
      taxa = taxa, degree = 0, within_links = 0,
      module_size = as.integer(table(module)[as.character(module)]),
      module_mean = 0, module_sd = 0,
      stringsAsFactors = FALSE, row.names = taxa
    )
    return(out)
  }

  ends <- igraph::ends(g, igraph::E(g), names = FALSE)
  from <- c(ends[, 1], ends[, 2])
  to   <- c(ends[, 2], ends[, 1])

  ## Conventional Zi-Pi uses edge counts; strength mode is optional.
  if (role_weights == "degree") {
    w <- rep(1, length(from))
  } else {
    ew <- NULL
    if ("weight" %in% igraph::edge_attr_names(g)) {
      ew0 <- suppressWarnings(as.numeric(igraph::E(g)$weight))
      if (length(ew0) == m && all(is.finite(ew0))) ew <- abs(ew0)
    }
    if (is.null(ew)) ew <- rep(1, m)
    w <- rep(ew, 2)
  }

  ## Accumulate each node's links/strength to every target module.
  ## key is unique for (node, target_module), so complexity is O(E).
  target_module <- module[to]
  key <- from + (target_module - 1L) * n

  agg <- rowsum(w, group = key, reorder = FALSE)
  agg_val <- as.numeric(agg[, 1])
  key_int <- as.integer(rownames(agg))
  node_id <- ((key_int - 1L) %% n) + 1L
  mod_id  <- ((key_int - 1L) %/% n) + 1L

  total <- numeric(n)
  tmp <- rowsum(agg_val, node_id, reorder = FALSE)
  total[as.integer(rownames(tmp))] <- tmp[, 1]

  sq <- numeric(n)
  tmp <- rowsum(agg_val^2, node_id, reorder = FALSE)
  sq[as.integer(rownames(tmp))] <- tmp[, 1]

  within <- numeric(n)
  own <- mod_id == module[node_id]
  if (any(own)) {
    tmp <- rowsum(agg_val[own], node_id[own], reorder = FALSE)
    within[as.integer(rownames(tmp))] <- tmp[, 1]
  }

  ## Pi = 1 - sum_s (k_is / k_i)^2
  p <- numeric(n)
  nz <- total > 0
  p[nz] <- 1 - sq[nz] / total[nz]^2
  p <- pmax(0, pmin(1, p))

  ## Zi = standardized within-module degree/strength.
  module_mean <- ave(within, module, FUN = mean)
  module_sd <- ave(within, module, FUN = stats::sd)
  z <- numeric(n)
  valid <- is.finite(module_sd) & module_sd > 0
  z[valid] <- (within[valid] - module_mean[valid]) / module_sd[valid]
  z[!is.finite(z)] <- 0

  role4 <- ifelse(
    z < zi_cut,
    ifelse(p < pi_cut, "Peripherals", "Connectors"),
    ifelse(p < pi_cut, "Module hubs", "Network hubs")
  )

  ## Original Guimera-Amaral seven-role partition.
  role7 <- ifelse(
    z < zi_cut,
    ifelse(p < .05, "ultra peripheral",
           ifelse(p < .62, "peripheral",
                  ifelse(p < .80, "non hub connector", "non hub kinless"))),
    ifelse(p < .30, "provincial hub",
           ifelse(p < .75, "connector hub", "kinless hub"))
  )

  out <- data.frame(
    z = z,
    module = module,
    p = p,
    roles = role4,
    role_7 = role7,
    taxa = taxa,
    degree = as.numeric(igraph::degree(g)),
    within_links = within,
    module_size = as.integer(table(module)[as.character(module)]),
    module_mean = module_mean,
    module_sd = module_sd,
    stringsAsFactors = FALSE,
    row.names = taxa
  )

  out
}


#' Multi-group Zi-Pi plot
#'
#' @param x Combined Zi-Pi data containing p, z and group.
#' @param zi_cut Zi hub threshold.
#' @param pi_cut Pi connector threshold.
#' @param label_key Label key-role nodes.
#' @export
facet.zipi.fast <- function(x, zi_cut = 2.5, pi_cut = 0.62, label_key = TRUE) {
  if (is.null(x) || !nrow(x)) return(NULL)
  if (!all(c("p","z","group") %in% names(x)))
    stop("x must contain p, z and group.")

  zones <- data.frame(
    xmin=c(0,pi_cut,0,pi_cut), xmax=c(pi_cut,1,pi_cut,1),
    ymin=c(-Inf,-Inf,zi_cut,zi_cut), ymax=c(zi_cut,zi_cut,Inf,Inf),
    role=c("Peripherals","Connectors","Module hubs","Network hubs")
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data=zones,
      ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=role),
      alpha=.16, inherit.aes=FALSE
    ) +
    ggplot2::geom_vline(xintercept=pi_cut,linetype=2) +
    ggplot2::geom_hline(yintercept=zi_cut,linetype=2) +
    ggplot2::geom_point(
      data=x,
      ggplot2::aes(x=p,y=z,color=factor(module)),
      size=2
    ) +
    ggplot2::facet_grid(.~group, scales="free") +
    ggplot2::guides(color="none") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x="Participation coefficient (Pi)",
      y="Within-module connectivity z-score (Zi)",
      fill="Topological roles"
    )

  if (label_key && requireNamespace("ggrepel", quietly=TRUE)) {
    lab <- x[x$z >= zi_cut | x$p >= pi_cut,,drop=FALSE]
    if (nrow(lab)) {
      if (!"label" %in% names(lab)) lab$label <- rownames(lab)
      p <- p + ggrepel::geom_text_repel(
        data=lab,
        ggplot2::aes(x=p,y=z,label=label,color=factor(module)),
        size=3, show.legend=FALSE
      )
    }
  }
  p
}



#The participation coefficient of a node measures how well a  node is distributed
# in the entire network. It is close to 1 if its links are uniformly
#distributed among all the modules and 0 if all its links are within its own module.

participation_coeffiecient <- function(mod.degree, total.degree){

  p <- NULL

  for(i in total.degree$taxa){

    ki <- subset(total.degree$total_links, total.degree$taxa==i)

    taxa.mod.degree <- subset(mod.degree$mod_links, mod.degree$taxa==i)

    p[i] <- 1 - (sum((taxa.mod.degree)**2)/ki**2)

  }

  p <- as.data.frame(p)

  return(p)

}


assign_module_roles <- function(zp){

  zp <- na.omit(zp)

  zp$roles <- rep(0, dim(zp)[1])

  outdf <- NULL

  for(i in 1:dim(zp)[1]){

    df <- zp[i, ]

    if(df$z < 2.5){ #non hubs

      if(df$p < 0.05){

        df$roles <- "ultra peripheral"

      }
      else if(df$p < 0.620){

        df$roles <- "peripheral"

      }
      else if(df$p < 0.80){

        df$roles <- "non hub connector"

      }
      else{

        df$roles <- "non hub kinless"

      }

    }
    else { # module hubs

      if(df$p < 0.3){

        df$roles <- "provincial hub"

      }
      else if(df$p < 0.75){

        df$roles <- "connector hub"

      }
      else {

        df$roles <- "kinless hub"

      }

    }

    if(is.null(outdf)){outdf <- df}else{outdf <- rbind(outdf, df)}

  }

  return(outdf)

}


plot_roles <- function(node.roles, roles.colors=NULL){

  x1<- c(0, 0.05, 0.62, 0.8, 0, 0.30, 0.75)
  x2<- c(0.05, 0.62, 0.80, 1,  0.30, 0.75, 1)
  y1<- c(-Inf,-Inf, -Inf, -Inf,  2.5, 2.5, 2.5)
  y2 <- c(2.5,2.5, 2.5, 2.5, Inf, Inf, Inf)

  lab <- c("ultra peripheral","peripheral" ,"non-hub connector","non-hub kinless","provincial"," hub connector","hub kinless")

  if(is.null(roles.colors)){roles.colors <- c("#E6E6FA", "#DCDCDC", "#F5FFFA", "#FAEBD7", "#EEE8AA", "#E0FFFF", "#F5F5DC")}

  p <- ggplot() + geom_rect(data=NULL, mapping=aes(xmin=x1, xmax=x2, ymin=y1,ymax=y2, fill=lab))

  p <- p + guides(fill=guide_legend(title="Topological roles"))

  p  <- p + scale_fill_manual(values = roles.colors)

  p <- p + geom_point(data=node.roles, aes(x=p, y=z,color=module)) + theme_bw()

  p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("Participation Coefficient")+ylab(" Within-module connectivity z-score")

  return(p)
}


# plot_roles3 = function(node.roles, roles.colors=NULL){
#
#   x1<- c(0, 0.62,0,0.62)
#   x2<- c( 0.62,1,0.62,1)
#   y1<- c(-Inf,2.5,2.5,-Inf)
#   y2 <- c(2.5,Inf,Inf,2.5)
#   #
#   lab <- c("peripheral",'Connectors','Module hubs','Network hubs')
#   #
#   if(is.null(roles.colors)){roles.colors <- c("#E6E6FA", "#DCDCDC","#F5FFFA", "#FAEBD7")}
#
#   p <- ggplot() + geom_rect(data=NULL,
#                             mapping=aes(xmin=x1, xmax=x2,ymin=y1,ymax=y2,fill = lab))
#   p
#   p <- p + guides(fill=guide_legend(title="Topological roles"))
#
#   p  <- p + scale_fill_manual(values = roles.colors)
#
#   p <- p + geom_point(data=node.roles, aes(x=p, y=z,color=module)) + theme_bw()
#
#   p<-p+theme(strip.background = element_rect(fill = "white"))+
#     xlab("Participation Coefficient")+ylab(" Within-module connectivity z-score")
#
#   return(p)
# }

plot_roles2 = function(node.roles, roles.colors=NULL){
  x1<- c(0, 0.62,0,0.62)
  x2<- c( 0.62,1,0.62,1)
  y1<- c(-Inf,2.5,2.5,-Inf)
  y2 <- c(2.5,Inf,Inf,2.5)
  lab <- c("peripheral",'Network hubs','Module hubs','Connectors')

  if(is.null(roles.colors)){ roles.colors <- c("#E6E6FA","#DCDCDC","#F5FFFA", "#FAEBD7")}

  p <- ggplot() + geom_rect(data=NULL,
                            mapping=aes(xmin=x1, xmax=x2,ymin=y1,ymax=y2,fill = lab))
  p
  p <- p + guides(fill=guide_legend(title="Topological roles"))

  p  <- p + scale_fill_manual(values = roles.colors)
  p <- p + geom_point(data=node.roles,aes(x=p, y=z,color=module)) + theme_bw()+
    guides(color= F)

  # 是否需要模块
  #p <- p + geom_point(data=node.roles,aes(x=p, y=z,color=module)) + theme_bw()

  #p <- p + geom_point(data=node.roles, aes(x=p, y=z)) + theme_bw()
  p<-p+theme(strip.background = element_rect(fill = "white"))+
    xlab("Participation Coefficient")+ylab(" Within-module connectivity z-score")
  return(p)
}
