# =================================================================================================
# DotPlot_allow_dups()
# -------------------------------------------------------------------------------------------------
# Purpose:
#   A wrapper around Seurat::DotPlot() that accepts duplicate gene names in the feature list.
#   Standard DotPlot() silently drops duplicates; this function re-introduces them so that the
#   same gene can appear multiple times (e.g. once per cluster group in a named list).
#
# What the function does:
#   • Accepts features as a plain character vector OR a named list of character vectors
#     (named list triggers Seurat-style faceting by group, identical to DotPlot behavior).
#   • Internally deduplicates the gene list before calling DotPlot(), then expands the
#     returned plot data to restore duplicated entries with unique internal suffixes.
#   • Axis labels display the original gene symbol; internal suffixes are hidden by default.
#
# Inputs:
#   object           : A Seurat object.
#   features         : Character vector or named list of character vectors. Duplicates allowed.
#   assay            : Assay to use (default: active assay).
#   cols             : Color scale passed to DotPlot (default: c("lightgrey", "blue")).
#   col.min / col.max: Color scale limits (default: -2.5 / 2.5).
#   dot.min          : Minimum dot size (default: 0).
#   dot.scale        : Dot scaling factor (default: 6).
#   idents           : Identity classes to include (default: all).
#   group.by         : Metadata column to group cells by.
#   split.by         : Metadata column to split dots by.
#   cluster.idents   : Whether to cluster identity classes (default: FALSE).
#   scale            : Whether to scale expression values (default: TRUE).
#   scale.by         : Scale dot size by "radius" or "size" (default: "radius").
#   scale.min / scale.max : Minimum/maximum scaled values.
#   dup_sep          : Internal string used to make duplicate gene IDs unique.
#                      Should not appear in any gene name (default: "___dup___").
#   show_dup_suffix  : If TRUE, display internal suffixed IDs on axis (default: FALSE).
#
# Outputs:
#   A ggplot2 object (same class as Seurat::DotPlot output), with duplicated genes
#   re-introduced in the plot data and displayed on the x-axis.
#
# Dependencies:
#   Seurat, ggplot2, grid (for facet spacing)
#
# Notes:
#   • Expression values and percent-expressed values for duplicated genes are identical
#     to the original entry — this is intentional; the duplication is purely visual.
#   • Faceting behavior (when features is a named list) mirrors Seurat::DotPlot exactly.
# =================================================================================================

DotPlot_allow_dups <- function(
    object,
    features,
    assay = NULL,
    cols = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA,
    dup_sep = "___dup___",   # internal suffix (hidden from axis labels)
    show_dup_suffix = FALSE  # set TRUE if you want to see OTX1___dup___2 etc.
) {
  stopifnot(inherits(object, "Seurat"))
  
  # ---- flatten + preserve grouping (Seurat-style facet behavior) ----
  feature_groups <- NULL
  if (is.list(features) || any(!is.na(names(features)))) {
    if (is.null(names(features))) names(features) <- paste0("Group", seq_along(features))
    feature_groups <- unlist(lapply(seq_along(features), function(i) {
      rep(names(features)[i], length(features[[i]]))
    }), use.names = FALSE)
    gene_vec <- unlist(features, use.names = FALSE)
  } else {
    gene_vec <- as.character(features)
  }
  
  # ---- make unique internal feature ids for duplicates ----
  # (keeps original order, just appends dup suffix where needed)
  feat_unique <- make.unique(gene_vec, sep = dup_sep)
  
  map_df <- data.frame(
    gene = gene_vec,
    feat_unique = feat_unique,
    stringsAsFactors = FALSE
  )
  if (!is.null(feature_groups)) {
    map_df$feature.groups <- feature_groups
  }
  
  # compute once per gene (Seurat default look)
  genes_unique <- unique(gene_vec)
  
  p <- Seurat::DotPlot(
    object = object,
    features = genes_unique,
    assay = assay,
    cols = cols,
    col.min = col.min,
    col.max = col.max,
    dot.min = dot.min,
    dot.scale = dot.scale,
    idents = idents,
    group.by = group.by,
    split.by = split.by,
    cluster.idents = cluster.idents,
    scale = scale,
    scale.by = scale.by,
    scale.min = scale.min,
    scale.max = scale.max
  ) + 
    theme_classic() + 
    theme(
      axis.title = element_blank(), 
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  
  # ---- expand plot data to re-introduce duplicates ----
  df <- p$data
  df$features.plot <- as.character(df$features.plot)
  
  df2 <- merge(
    x = df,
    y = map_df,
    by.x = "features.plot",
    by.y = "gene",
    all.x = TRUE,
    sort = FALSE
  )
  
  # set x as the unique internal ids (avoids duplicated factor levels)
  df2$features.plot <- factor(df2$feat_unique, levels = map_df$feat_unique)
  
  # add Seurat-like facetting by groups if user provided a list
  if (!is.null(feature_groups)) {
    df2$feature.groups <- factor(df2$feature.groups, levels = unique(feature_groups))
  }
  
  p$data <- df2
  
  # show original gene symbols on the axis (hide internal suffix)
  if (!show_dup_suffix) {
    p <- p + ggplot2::scale_x_discrete(labels = function(x) sub(paste0(dup_sep, "\\d+$"), "", x))
  }
  
  # facet like Seurat does, with free_x + free space (variable widths by #genes)
  if (!is.null(feature_groups)) {
    p <- p +
      ggplot2::facet_grid(
        facets = ~feature.groups,
        scales = "free_x",
        space  = "free_x",
        switch = "y"
      ) +
      ggplot2::theme(
        panel.spacing = grid::unit(1, "lines"),
        strip.background = ggplot2::element_blank()
      )
  }
  
  return(p)
}
