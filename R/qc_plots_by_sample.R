# =================================================================================================
# qc_plots_by_sample()
# -------------------------------------------------------------------------------------------------
# Purpose:
#   Automated QC thresholding and visualization for single-cell datasets on a per-sample basis.
#   For each sample, this function:
#     • Computes data-driven QC minimum thresholds using density valley detection
#       (with optional Otsu/GMM assistance and guardrails).
#     • Applies fixed QC thresholds for metrics that do not use valley-based minima.
#     • Generates a one-page PDF containing:
#           - density plots (with detected valleys)
#           - violin plots with threshold lines
#           - scatter plots (log2 metrics vs. QC metrics, and pairwise log2 metrics)
#     • Records all thresholds and chosen strategies in a summary CSV file.
#
# Key Features:
#   • Fully sample-aware: thresholds are computed independently for each sample.
#   • Valley detection methods:
#         - Identifies candidate valleys in density(log2(metric)).
#         - Score-based ranking combines valley depth, effect size, and class balance.
#         - Incorporates Otsu threshold and optional GMM mixture boundaries.
#   • Guardrails:
#         - Rejects minima that are too high (>90th percentile or drop too many cells).
#         - Falls back to low quantile threshold when no valley is suitable.
#   • Flexible rule system:
#         - Global default rules, global per-metric rules, or per-sample overrides.
#         - Strategies supported: auto, first, deepest, nth, closest_to, left_of, right_of.
#         - Manual raw or log2 minima override all automated rules.
#   • Optional per-metric min/max clamps (floors_raw, caps_raw).
#   • Produces reproducible, publication-ready QC summaries.
#
# Inputs:
#   seurat_obj            : A Seurat object with QC metrics in meta.data.
#   sample_col            : Column name defining sample identity (e.g. orig.ident).
#   log2_qc_cols          : Metrics processed on log2 scale with valley detection (e.g., nCount_RNA).
#   qc_cols               : Metrics using fixed thresholds only (e.g., percent.mt).
#   path                  : Output directory for PDFs and summary CSV.
#   group.by              : Optional meta.data column for coloring scatter plots.
#   upper_quantile        : Raw upper bound for max threshold (e.g. 99.5%).
#   lower_quantile_no_valley :
#                           Raw min threshold when no valley is found (e.g. 0.5%).
#   max_min_quantile      : Guardrail: reject minima >= this quantile.
#   max_drop_fraction     : Guardrail: reject minima that remove too many cells.
#   mt_max                : Default percent.mt maximum if none is provided.
#   fixed_thresholds      : Named list of min/max values for fixed QC metrics.
#   valley_rules          : Global rules for valley selection per metric.
#   valley_rules_by_sample:
#                           Per-sample overrides for metric-specific valley rules.
#   floors_raw / caps_raw : Optional metric-specific min/max clamps for valley minima.
#   pdf_width, pdf_height : PDF size parameters.
#   plots_ncol            : Number of columns per PDF page.
#
# Outputs:
#   • A set of single-page PDFs (one per sample) containing all QC plots.
#   • "QC_thresholds_summary.csv" summarizing:
#         - sample
#         - metric
#         - min / max thresholds
#         - chosen valley (log2 value)
#         - method used (strategy / fallback)
#   • Invisibly returns a list: 
#         $thresholds (data.frame), $pdfs (file paths)
#
# Dependencies:
#   Seurat, dplyr, ggplot2, purrr, tibble, readr, patchwork
#   Optional: mclust (for Gaussian mixture valley suggestions)
#
# Notes:
#   • Designed for robust automated QC in large multi-sample scRNA-seq/scATAC-like datasets.
#   • Highly configurable — supports per-sample customization, manual minima, and safe fallback logic.
# =================================================================================================


# ---- QC plotter with per-sample forced valleys, manual mins, guards, one-page layout ----
# deps: Seurat, dplyr, ggplot2, purrr, tibble, readr, patchwork  (optional: mclust for GMM)
qc_plots_by_sample <- function(
    seurat_obj,
    sample_col,                                # e.g. "orig.ident" / "donor"
    log2_qc_cols = c("nCount_RNA","nFeature_RNA"),
    qc_cols      = c("percent.mt"),
    path         = "qc_plots",
    group.by     = NULL,                       # e.g. "CellType_AllCells" (color in scatters)
    upper_quantile = 0.995,                    # RAW scale upper bound (e.g. 99.5%)
    lower_quantile_no_valley = 0.005,          # RAW scale min when no valleys (0.5%)
    max_min_quantile = 0.90,                   # too-high rule 1: cap by percentile
    max_drop_fraction = 0.70,                  # too-high rule 2: cap by % dropped
    mt_max = 5,                                # convenience default for percent.mt
    fixed_thresholds = list(),                 # named list: metric -> list(min=?, max=?)
    valley_rules = list(.default = list(strategy = "auto")),     # global: strategies
    valley_rules_by_sample = list(),           # per-sample overrides: sample -> list(metric/.default -> rule)
    floors_raw = list(),                       # optional min clamps (RAW)
    caps_raw   = list(),                       # optional max clamps (RAW)
    pdf_width = 11, pdf_height = 8.5,          # one page per sample
    plots_ncol = 2                             # number of columns on the page
){
  stopifnot(sample_col %in% colnames(seurat_obj@meta.data))
  # Requires: Seurat, dplyr, ggplot2, purrr, tibble, readr, patchwork
  # (optional: mclust for Gaussian mixture valley suggestions)
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  .to_log2x <- function(raw) log2(raw + 1)
  .to_raw   <- function(log2x) (2^log2x) - 1
  
  # ---------- valley helpers ----------
  .valley_indices <- function(y) which(diff(sign(diff(y))) == 2) + 1L
  .valley_prominence <- function(x, y){
    v <- .valley_indices(y); if(!length(v)) return(numeric(0))
    prom <- numeric(length(v))
    for (i in seq_along(v)){
      lp <- { left <- y[1:v[i]]; idx <- which(diff(sign(diff(left))) == -2) + 1L
      if (length(idx)) idx[length(idx)] else 1L }
      rp <- { right <- y[v[i]:length(y)]; idx <- which(diff(sign(diff(right))) == -2) + 1L
      if (length(idx)) (idx[1] + v[i] - 1L) else length(y) }
      prom[i] <- min(y[lp], y[rp]) - y[v[i]]
    }
    prom
  }
  
  # ---------- Otsu / (optional) GMM helpers ----------
  .otsu_threshold <- function(x, nbins = 256){
    x <- x[is.finite(x)]
    h <- hist(x, breaks = nbins, plot = FALSE)
    p <- h$counts / sum(h$counts)
    om <- cumsum(p)
    mu <- cumsum(p * h$mids)
    mu_t <- sum(p * h$mids)
    sb2 <- (mu_t*om - mu)^2 / pmax(om*(1-om), 1e-12)
    h$mids[ which.max(sb2) ]
  }
  .gmm_boundaries <- function(x){
    if (!requireNamespace("mclust", quietly = TRUE)) return(numeric(0))
    fit <- try(mclust::Mclust(x, G = 2:4, verbose = FALSE), silent = TRUE)
    if (inherits(fit, "try-error") || is.null(fit$parameters$mean)) return(numeric(0))
    mu  <- as.numeric(fit$parameters$mean)
    sig <- sqrt(as.numeric(fit$parameters$variance$sigmasq))
    pi_ <- as.numeric(fit$parameters$pro)
    ord <- order(mu); mu <- mu[ord]; sig <- sig[ord]; pi_ <- pi_[ord]
    bounds <- c()
    for (i in seq_len(length(mu)-1)){
      f <- function(z) pi_[i]*dnorm(z, mu[i], sig[i]) - pi_[i+1]*dnorm(z, mu[i+1], sig[i+1])
      lo <- min(mu[i], mu[i+1]); hi <- max(mu[i], mu[i+1])
      if (sign(f(lo)) == sign(f(hi))) next
      r <- try(uniroot(f, c(lo, hi))$root, silent = TRUE)
      if (!inherits(r, "try-error")) bounds <- c(bounds, r)
    }
    bounds
  }
  
  .score_split <- function(x_log2, t_log2, dens_x, dens_y){
    if (!is.finite(t_log2)) return(-Inf)
    y_at_t <- approx(dens_x, dens_y, xout = t_log2, rule = 2)$y
    s1 <- 1 - (y_at_t / max(dens_y))
    left <- x_log2 <= t_log2; n <- length(x_log2)
    n1 <- sum(left); n2 <- n - n1
    if (n1 < 5 || n2 < 5) return(-Inf)
    m1 <- mean(x_log2[left]);  m2 <- mean(x_log2[!left])
    sd1 <- sd(x_log2[left]);   sd2 <- sd(x_log2[!left])
    sp  <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/pmax(n1+n2-2,1))
    d   <- if (is.finite(sp) && sp>0) abs(m2 - m1)/sp else 0
    s2  <- d/(d+1)
    bal <- 2*min(n1/n, n2/n)
    0.5*s1 + 0.35*s2 + 0.15*bal
  }
  
  .candidate_valleys_ranked <- function(vec_raw){
    x_log2 <- .to_log2x(vec_raw)
    d <- density(x_log2, na.rm = TRUE)
    v_idx <- .valley_indices(d$y)
    v_x   <- if (length(v_idx)) d$x[v_idx] else numeric(0)
    cand_x <- v_x
    t_otsu <- .otsu_threshold(x_log2)
    if (length(v_x)) cand_x <- unique(c(cand_x, v_x[ which.min(abs(v_x - t_otsu)) ]))
    t_gmm <- .gmm_boundaries(x_log2)
    if (length(t_gmm) && length(v_x)) {
      nearest <- v_x[ apply(abs(outer(v_x, t_gmm, "-")), 1, which.min) ]
      cand_x <- unique(c(cand_x, nearest))
    }
    scores <- if (length(cand_x)) vapply(cand_x, function(cx) .score_split(x_log2, cx, d$x, d$y), numeric(1)) else numeric(0)
    ord <- if (length(scores)) order(scores, decreasing = TRUE) else integer(0)
    list(d = d, v_idx = v_idx, all_valley_x = v_x, cand_x = if(length(ord)) cand_x[ord] else cand_x)
  }
  
  .choose_valley_rule <- function(d, rule){
    v_idx <- .valley_indices(d$y); if (!length(v_idx)) return(NA_integer_)
    strat <- (rule$strategy %||% "deepest")
    if (strat == "first") return(v_idx[1])
    if (strat == "deepest"){
      prom <- .valley_prominence(d$x, d$y); if (!length(prom)) return(v_idx[1])
      return(v_idx[ which.max(prom) ])
    }
    if (strat == "nth"){
      n <- as.integer(rule$n %||% 1L); n <- max(1L, n)
      if (n > length(v_idx)) return(v_idx[length(v_idx)]) else return(v_idx[n])
    }
    if (strat == "closest_to"){
      tx <- .to_log2x(as.numeric(rule$value)); return(v_idx[ which.min(abs(d$x[v_idx] - tx)) ])
    }
    if (strat == "right_of"){
      tx <- .to_log2x(as.numeric(rule$value)); cand <- v_idx[ d$x[v_idx] >= tx ]; if (length(cand)) cand[1] else v_idx[length(v_idx)]
    } else if (strat == "left_of"){
      tx <- .to_log2x(as.numeric(rule$value)); cand <- v_idx[ d$x[v_idx] <= tx ]; if (length(cand)) cand[length(cand)] else v_idx[1]
    } else v_idx[1]
  }
  
  .is_too_high <- function(min_raw, vec_raw){
    if (!is.finite(min_raw)) return(TRUE)
    qcap <- as.numeric(stats::quantile(vec_raw, probs = max_min_quantile, na.rm = TRUE))
    drop_frac <- mean(vec_raw < min_raw, na.rm = TRUE)
    (min_raw >= qcap) || (drop_frac > max_drop_fraction)
  }
  
  .compute_thresholds <- function(vec_raw, metric, rule_for_metric){
    vec_raw <- vec_raw[is.finite(vec_raw)]
    
    # --- manual minima (bypass valley picking if provided) ---
    if (length(vec_raw) && (!is.null(rule_for_metric$raw_min) || !is.null(rule_for_metric$log2_min))){
      if (!is.null(rule_for_metric$raw_min)) {
        min_thr <- as.numeric(rule_for_metric$raw_min)
        v_x     <- .to_log2x(min_thr)
        picked_from <- "manual_raw_min"
      } else {
        v_x     <- as.numeric(rule_for_metric$log2_min)
        min_thr <- .to_raw(v_x)
        picked_from <- "manual_log2_min"
      }
      max_thr <- stats::quantile(vec_raw, probs = upper_quantile, na.rm = TRUE)
      fl <- as.numeric(floors_raw[[metric]] %||% NA_real_)
      cp <- as.numeric(caps_raw[[metric]]   %||% NA_real_)
      if (is.finite(fl)) min_thr <- max(min_thr, fl, na.rm = TRUE)
      if (is.finite(cp)) min_thr <- min(min_thr, cp, na.rm = TRUE)
      return(list(min = as.numeric(min_thr), max = as.numeric(max_thr), valley_x = as.numeric(v_x),
                  all_valleys = integer(0), all_valley_x = numeric(0),
                  dens_x = density(.to_log2x(vec_raw))$x, dens_y = density(.to_log2x(vec_raw))$y,
                  picked_from = picked_from))
    }
    
    if (!length(vec_raw)) {
      return(list(min = NA_real_, max = NA_real_, valley_x = NA_real_,
                  all_valleys = integer(0), all_valley_x = numeric(0),
                  dens_x = numeric(0), dens_y = numeric(0),
                  picked_from = "none"))
    }
    
    ranked <- .candidate_valleys_ranked(vec_raw)
    d <- ranked$d; v_idx <- ranked$v_idx
    
    # no valleys -> 0.5% fallback
    if (!length(v_idx)){
      min_thr <- as.numeric(stats::quantile(vec_raw, probs = lower_quantile_no_valley, na.rm = TRUE))
      max_thr <- as.numeric(stats::quantile(vec_raw, probs = upper_quantile,           na.rm = TRUE))
      v_x     <- .to_log2x(min_thr)
      fl <- as.numeric(floors_raw[[metric]] %||% NA_real_)
      cp <- as.numeric(caps_raw[[metric]]   %||% NA_real_)
      if (is.finite(fl)) min_thr <- max(min_thr, fl, na.rm = TRUE)
      if (is.finite(cp)) min_thr <- min(min_thr, cp, na.rm = TRUE)
      return(list(min = min_thr, max = max_thr, valley_x = v_x,
                  all_valleys = integer(0), all_valley_x = numeric(0),
                  dens_x = d$x, dens_y = d$y, picked_from = "fallback_q0.5%"))
    }
    
    # resolve strategy & allow_adjust
    allow_adjust <- isTRUE(rule_for_metric$allow_adjust %||% TRUE)
    strat <- (rule_for_metric$strategy %||% "auto")
    
    if (strat == "auto"){
      v_x <- if (length(ranked$cand_x)) ranked$cand_x[1] else NA_real_
    } else {
      picked <- .choose_valley_rule(d, rule_for_metric)
      v_x <- if (is.na(picked)) NA_real_ else d$x[picked]
    }
    
    min_thr <- if (is.na(v_x)) NA_real_ else .to_raw(v_x)
    picked_from <- strat
    
    # adjust only if allowed
    if (allow_adjust){
      if (.is_too_high(min_thr, vec_raw)){
        ok <- FALSE
        if (length(ranked$cand_x)){
          for (cx in ranked$cand_x){
            test_min <- .to_raw(cx)
            if (!.is_too_high(test_min, vec_raw)){
              v_x <- cx; min_thr <- test_min; picked_from <- paste0(strat,"_adj")
              ok <- TRUE; break
            }
          }
        }
        if (!ok){
          min_thr <- as.numeric(stats::quantile(vec_raw, probs = lower_quantile_no_valley, na.rm = TRUE))
          v_x <- .to_log2x(min_thr)
          picked_from <- "fallback_q0.5%"
        }
      }
    } else {
      picked_from <- paste0(strat, "_forced")
    }
    
    max_thr <- stats::quantile(vec_raw, probs = upper_quantile, na.rm = TRUE)
    fl <- as.numeric(floors_raw[[metric]] %||% NA_real_)
    cp <- as.numeric(caps_raw[[metric]]   %||% NA_real_)
    if (is.finite(fl)) min_thr <- max(min_thr, fl, na.rm = TRUE)
    if (is.finite(cp)) min_thr <- min(min_thr, cp, na.rm = TRUE)
    
    list(min = as.numeric(min_thr),
         max = as.numeric(max_thr),
         valley_x = as.numeric(v_x),
         all_valleys = v_idx,
         all_valley_x = d$x[v_idx],
         dens_x = d$x, dens_y = d$y,
         picked_from = picked_from)
  }
  
  # ---------- plotting ----------
  .density_plot <- function(metric, thr_info){
    df <- tibble(x = thr_info$dens_x, y = thr_info$dens_y)
    p <- ggplot(df, aes(x, y)) +
      geom_line() +
      labs(title = paste0(metric, " (density on log2)"),
           x = "log2(value + 1)", y = "Density") +
      theme_classic(base_size = 8)
    if (length(thr_info$all_valleys)){
      vtab <- tibble(ix = thr_info$all_valleys,
                     vx = thr_info$all_valley_x,
                     vy = df$y[thr_info$all_valleys],
                     lab = seq_along(thr_info$all_valleys))
      p <- p + geom_point(data = vtab, aes(vx, vy), size = 2, alpha = 0.8) +
        geom_text(data = vtab, aes(vx, vy, label = lab), vjust = -0.7, size = 3)
    }
    if (is.finite(thr_info$valley_x)){
      vy <- df$y[ which.min(abs(df$x - thr_info$valley_x)) ]
      p <- p + geom_point(aes(x = thr_info$valley_x, y = vy), color = "red", size = 2.5) +
        geom_vline(xintercept = thr_info$valley_x, linetype = "dashed", color = "red")
    }
    p
  }
  
  # ---------- main loop ----------
  fixed_thr_eff <- fixed_thresholds
  if (is.null(fixed_thr_eff$`percent.mt`) && "percent.mt" %in% qc_cols) {
    fixed_thr_eff$`percent.mt` <- list(max = mt_max)
  }
  
  samples <- sort(unique(seurat_obj@meta.data[[sample_col]]))
  out_rows <- list(); pdf_files <- character(0)
  
  for (s in samples){
    cells_s <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[sample_col]] == s]
    if (!length(cells_s)) next
    sobj_s <- subset(seurat_obj, cells = cells_s)
    meta_s <- sobj_s@meta.data
    
    thr_this <- list()
    plots <- list()
    
    .resolve_rule <- function(metric){
      rs <- valley_rules_by_sample[[as.character(s)]]
      rs_metric  <- if (!is.null(rs)) rs[[metric]] else NULL
      rs_default <- if (!is.null(rs)) rs[[".default"]] else NULL
      global_metric  <- valley_rules[[metric]]
      global_default <- valley_rules[[".default"]]
      rs_metric %||% rs_default %||% global_metric %||% global_default %||% list(strategy="auto")
    }
    
    # A) thresholds for log2_qc_cols
    for (col in log2_qc_cols){
      if (!col %in% colnames(meta_s)) next
      rule <- .resolve_rule(col)
      ti <- .compute_thresholds(meta_s[[col]], metric = col, rule_for_metric = rule)
      thr_this[[col]] <- list(min = ti$min, max = ti$max, valley_x = ti$valley_x,
                              dens = ti, type = paste0(ti$picked_from, "+", upper_quantile*100, "pctl"))
      plots[[length(plots)+1]] <- .density_plot(col, ti)
    }
    
    # B) fixed thresholds for qc_cols (generic)
    for (col in qc_cols){
      if (!col %in% colnames(meta_s)) next
      fx <- fixed_thr_eff[[col]]
      tmin <- tmax <- NA_real_
      if (is.numeric(fx)) { tmax <- as.numeric(fx) }
      if (is.list(fx))   { tmin <- as.numeric(fx$min %||% NA_real_); tmax <- as.numeric(fx$max %||% NA_real_) }
      thr_this[[col]] <- list(min = tmin, max = tmax, valley_x = NA_real_, dens = NULL, type = "fixed")
    }
    
    # C) violins with threshold lines
    add_thr_lines <- function(p, col){
      tmin <- thr_this[[col]]$min; tmax <- thr_this[[col]]$max
      if (is.finite(tmin)) p <- p + geom_hline(yintercept = tmin, linetype = "dashed", color = "red")
      if (is.finite(tmax)) p <- p + geom_hline(yintercept = tmax, linetype = "dashed", color = "red")
      p
    }
    for (col in unique(c(log2_qc_cols, qc_cols))){
      if (!col %in% colnames(meta_s)) next
      p <- VlnPlot(sobj_s, features = col, pt.size = 0) +
        theme_classic(base_size = 8) + theme(legend.position = "none", axis.title.x = element_blank()) +
        labs(title = paste0(col))
      plots[[length(plots)+1]] <- add_thr_lines(p, col)
    }
    
    # D) scatters
    do_FS <- function(f1, f2){
      if (!all(c(f1,f2) %in% colnames(meta_s))) return(NULL)
      p <- FeatureScatter(sobj_s, feature1 = f1, feature2 = f2, group.by = group.by) +
        labs(title = paste0(f1, " vs ", f2)) + theme_classic(base_size = 8)
      t1 <- thr_this[[f1]]; t2 <- thr_this[[f2]]
      if (!is.null(t1)) {
        if (is.finite(t1$min)) p <- p + geom_vline(xintercept = t1$min, linetype = "dashed", color = "red")
        if (is.finite(t1$max)) p <- p + geom_vline(xintercept = t1$max, linetype = "dashed", color = "red")
      }
      if (!is.null(t2)) {
        if (is.finite(t2$min)) p <- p + geom_hline(yintercept = t2$min, linetype = "dashed", color = "red")
        if (is.finite(t2$max)) p <- p + geom_hline(yintercept = t2$max, linetype = "dashed", color = "red")
      }
      p
    }
    for (xname in log2_qc_cols) for (yname in qc_cols) {
      p <- do_FS(xname, yname); if (!is.null(p)) plots[[length(plots)+1]] <- p
    }
    if (length(log2_qc_cols) >= 2){
      combs <- combn(log2_qc_cols, 2, simplify = FALSE)
      for (cmb in combs) {
        p <- do_FS(cmb[[1]], cmb[[2]]); if (!is.null(p)) plots[[length(plots)+1]] <- p
      }
    }
    
    # ---- write single-page PDF ----
    pdf_file <- file.path(path, sprintf("QC_%s_%s.pdf", sample_col, as.character(s)))
    grDevices::pdf(pdf_file, width = pdf_width, height = pdf_height)
    if (length(plots) == 0) {
      grid::grid.newpage(); grid::grid.text(sprintf("%s: no plots", s))
    } else {
      page <- patchwork::wrap_plots(plots, ncol = plots_ncol)
      page <- page + patchwork::plot_annotation(
        title = sprintf("%s = %s", sample_col, as.character(s)),
        theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
      )
      print(page)
    }
    grDevices::dev.off()
    pdf_files <- c(pdf_files, pdf_file)
    
    # summary rows
    for (nm in names(thr_this)){
      out_rows[[length(out_rows)+1]] <- tibble(
        sample_col = sample_col, sample = as.character(s),
        metric = nm,
        min_thr = as.numeric(thr_this[[nm]]$min),
        max_thr = as.numeric(thr_this[[nm]]$max),
        valley_log2 = as.numeric(thr_this[[nm]]$valley_x),
        method = thr_this[[nm]]$type
      )
    }
  }
  
  summary_df <- dplyr::bind_rows(out_rows)
  out_csv <- file.path(path, "QC_thresholds_summary.csv")
  readr::write_csv(summary_df, out_csv)
  message("QC plots written: ", length(pdf_files))
  message("Threshold summary: ", out_csv)
  invisible(list(thresholds = summary_df, pdfs = pdf_files))
}
