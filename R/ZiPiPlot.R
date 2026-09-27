#' Fast and robust Zi-Pi topological roles
#'
#' Drop-in replacement for the previous ZiPiPlot/module.roles implementation.
#' Zi follows within-module degree z-score; Pi follows Guimera & Amaral (2005).
#' Computation is O(E) after community detection instead of O(V^2).
#'
#' @param igraph igraph object.
#' @param method Community method: cluster_fast_greedy, cluster_walktrap,
#'   cluster_edge_betweenness, cluster_spinglass, or cluster_louvain.
#' @param membership Optional precomputed module membership. If supplied,
#'   community detection is skipped.
#' @param module_weights Use absolute edge weights for community detection.
#' @param role_weights "degree" (recommended/default) or "strength".
#' @param zi_cut Zi threshold for hubs.
#' @param pi_cut Pi threshold for the simplified four-role scheme.
#' @param label_key Label connectors/module hubs/network hubs.
#' @return ZiPiPlot returns list(plot, data); module.roles returns data.frame.
#' @export
ZiPiPlot <- function(
    igraph,
    method = "cluster_fast_greedy",
    membership = NULL,
    module_weights = TRUE,
    role_weights = c("degree", "strength"),
    zi_cut = 2.5,
    pi_cut = 0.62,
    label_key = TRUE
) {
  role_weights <- match.arg(role_weights)
  g <- igraph

  if (!inherits(g, "igraph")) stop("igraph must be an igraph object.")
  if (igraph::vcount(g) == 0L) stop("igraph has no vertices.")

  if (is.null(igraph::V(g)$name))
    igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))

  ## Zi-Pi is defined here for undirected ecological networks.
  if (igraph::is_directed(g))
    g <- igraph::as_undirected(g, mode = "collapse")

  ## Remove loops/multiple edges; keep a numeric correlation/weight if available.
  if (igraph::any_loop(g) || any(igraph::count_multiple(g) > 1L)) {
    comb <- list(weight = function(x) {
      x <- suppressWarnings(as.numeric(x))
      if (!length(x) || all(!is.finite(x))) return(1)
      x <- x[is.finite(x)]
      x[which.max(abs(x))]
    }, "ignore")
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                          edge.attr.comb = comb)
  }

  n <- igraph::vcount(g)
  m <- igraph::ecount(g)

  ## Robust numeric edge weights.
  ew <- NULL
  if ("weight" %in% igraph::edge_attr_names(g) && m > 0L) {
    ew0 <- suppressWarnings(as.numeric(igraph::E(g)$weight))
    if (length(ew0) == m && all(is.finite(ew0))) ew <- ew0
  }
  if (is.null(ew) && m > 0L) ew <- rep(1, m)

  ## -------- module detection --------
  if (is.null(membership)) {
    if (m == 0L) {
      membership <- seq_len(n)
    } else {
      cw <- if (module_weights) abs(ew) else NULL
      if (!is.null(cw) && (!length(cw) || all(cw == 0))) cw <- NULL

      fc <- switch(
        method,
        cluster_fast_greedy =
          igraph::cluster_fast_greedy(g, weights = cw),
        cluster_walktrap =
          igraph::cluster_walktrap(g, weights = cw),
        cluster_louvain =
          igraph::cluster_louvain(g, weights = cw),
        cluster_edge_betweenness = {
          ## igraph interprets these weights as distances.
          dd <- if (is.null(cw)) NULL else 1 / pmax(cw, .Machine$double.eps)
          igraph::cluster_edge_betweenness(g, weights = dd)
        },
        cluster_spinglass =
          igraph::cluster_spinglass(g, weights = cw),
        stop("Unsupported community method: ", method)
      )
      membership <- igraph::membership(fc)
    }
  }

  if (length(membership) != n)
    stop("membership length must equal vcount(igraph).")

  ## Normalize arbitrary module labels to compact integer IDs.
  membership <- match(as.character(membership), unique(as.character(membership)))
  names(membership) <- igraph::V(g)$name
  igraph::V(g)$module <- as.character(membership)

  roles <- .zipi_roles_fast(
    g,
    membership = membership,
    role_weights = role_weights,
    zi_cut = zi_cut,
    pi_cut = pi_cut
  )

  roles$label <- if (label_key && nrow(roles))
    ifelse(roles$z >= zi_cut | roles$p >= pi_cut, roles$taxa, "") else ""

  ## Four-zone plot retained for backwards compatibility.
  zones <- data.frame(
    xmin = c(0, pi_cut, 0, pi_cut),
    xmax = c(pi_cut, 1, pi_cut, 1),
    ymin = c(-Inf, -Inf, zi_cut, zi_cut),
    ymax = c(zi_cut, zi_cut, Inf, Inf),
    role = c("Peripherals", "Connectors", "Module hubs", "Network hubs")
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = zones,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = role),
      alpha = 0.16, inherit.aes = FALSE
    ) +
    ggplot2::geom_vline(xintercept = pi_cut, linetype = 2) +
    ggplot2::geom_hline(yintercept = zi_cut, linetype = 2) +
    ggplot2::geom_point(
      data = roles,
      ggplot2::aes(x = p, y = z, color = factor(module)),
      size = 2
    ) +
    ggplot2::theme_bw() +
    ggplot2::guides(color = "none") +
    ggplot2::labs(
      x = "Participation coefficient (Pi)",
      y = "Within-module connectivity z-score (Zi)",
      fill = "Topological roles"
    )

  if (label_key && requireNamespace("ggrepel", quietly = TRUE) &&
      any(nzchar(roles$label))) {
    p <- p + ggrepel::geom_text_repel(
      data = roles[nzchar(roles$label), , drop = FALSE],
      ggplot2::aes(x = p, y = z, label = label, color = factor(module)),
      size = 3, show.legend = FALSE
    )
  }

  list(p, roles)
}

