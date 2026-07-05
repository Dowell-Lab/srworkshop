# =============================================================================
# ADVANCED track - Script 08: GSEA over GO Biological Process (gseGO / fgsea)
# -----------------------------------------------------------------------------
# WHAT THIS TEACHES
#   How to turn a per-gene differential-expression table into PATHWAY-level
#   biology using Gene Set Enrichment Analysis (GSEA). Instead of asking "is
#   this one gene significant?", GSEA asks "are the genes of a known pathway
#   collectively shifted up or down?" - which is far more robust to noise and
#   far easier to interpret.
#
#   Pipeline (mirrors the Spencer Lab DE+GSEA pipeline):
#     1. DE per cell type / condition with a ranked, UNFILTERED gene list
#        (no min.pct, no logfc cutoff - GSEA needs the whole ranking).
#     2. Rank statistic = avg_log2FC * -log10(p_val_adj)
#        (p_val_adj == 0 replaced by .Machine$double.xmin so -log10 stays finite).
#     3. Map SYMBOL -> ENTREZID (bitr / org.Hs.eg.db), then gseGO(ont="BP",
#        by="fgsea"). Return ALL terms (pvalueCutoff = 1) and filter downstream.
#     4. Visualize: two-panel NES dot plot (Activated vs Suppressed) + volcano.
#
# DATASET: GATA1 (GSE271399), annotated object from gata1/03. We compare the
#   two constructs (GATA1s vs wtGATA1) WITHIN each SingleR cell type.
#
# INPUT : OUT_DIR/gata1_combined_annotated.rds
# STYLE : ports the lab's helper functions (save_plot, pretty_label,
#         run_gsea_go_bp, make_volcano, plot_gsea_two_panel) with GATA1 metadata
#         (construct / SingleR labels) in place of the template's Treatment /
#         seurat_clusters, adapted to course conventions (read-only inputs,
#         write to OUT_DIR).
# =============================================================================

source("../00_paths_and_setup.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(forcats)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})
# presto makes Seurat DE run in minutes instead of hours (optional but great):
#   remotes::install_github("immunogenomics/presto")
suppressWarnings(suppressMessages(requireNamespace("presto", quietly = TRUE)))

# Keep memory sane on the cluster (as in the lab's original script):
library(future)
plan("sequential")
options(future.globals.maxSize = Inf)

# ---- Parameter block --------------------------------------------------------
io <- list(
  sobj_rds     = file.path(OUT_DIR, "gata1_combined_annotated.rds"),
  out_dir      = OUT_DIR,
  ident_col    = "SingleR_labels_other",  # cell type we test WITHIN
  contrast     = "construct",             # factor being compared
  test_level   = "GATA1s",                # numerator (ident.1)
  ref_level    = "wtGATA1",               # denominator (ident.2)
  assay        = "RNA",
  # DE parameters: pass the FULL, unfiltered ranking to GSEA
  de_test      = "wilcox",
  min_pct      = 0,
  logfc_thresh = 0,
  # GSEA parameters
  gsea_min     = 15,
  gsea_max     = 500,
  gsea_eps     = 1e-20,
  # Volcano thresholds
  volcano_lfc  = 1.5,
  volcano_padj = 0.05,
  # Only run GSEA for cell types with at least this many cells per side
  min_cells    = 30
)
outdir_gsea <- file.path(io$out_dir, paste0("advanced_gsea_", io$contrast))
dir.create(outdir_gsea, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# HELPER FUNCTIONS  (ported from the Spencer Lab DE+GSEA pipeline)
# =============================================================================

# ---- save_plot(): save a ggplot to PDF or PNG -------------------------------
save_plot <- function(plot, filename, path = ".", height = 6, width = 8,
                      resolution = 150, format = "png") {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  file_path <- file.path(path, paste0(filename, ".", format))
  if (format == "png") {
    png(file = file_path, height = height, width = width, units = "in", res = resolution)
    print(plot); dev.off()
  } else {
    pdf(file = file_path, height = height, width = width, useDingbats = FALSE)
    print(plot); dev.off()
  }
  message("Plot saved: ", file_path)
  invisible(file_path)
}

# ---- clean_name(): safe filename fragment -----------------------------------
clean_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

# ---- pretty_label(): strip GO/KEGG/etc. prefixes for plotting ---------------
pretty_label <- function(x, prefix = NULL) {
  if (!is.null(prefix)) {
    x <- sub(paste0("^", prefix, "_?"), "", x)
  } else {
    x <- sub(paste0("^(HALLMARK|GOBP|GOCC|GOMF|KEGG|REACTOME|BIOCARTA|WP|PID)_"),
             "", x)
  }
  gsub("_", " ", x)
}

# ---- to_df(): GSEA result object -> tibble ----------------------------------
to_df <- function(gsea_obj) {
  if (is.null(gsea_obj) || is.null(gsea_obj@result)) return(tibble::tibble())
  tibble::as_tibble(gsea_obj@result)
}

# ---- run_gsea_go_bp(): the core GSEA function -------------------------------
# Rank statistic: avg_log2FC * -log10(p_val_adj); p_val_adj == 0 -> double.xmin
run_gsea_go_bp <- function(de_df,
                           minGSSize     = 15,
                           maxGSSize     = 500,
                           eps           = 1e-20,
                           pAdjustMethod = "BH",
                           seed          = TRUE,
                           verbose       = FALSE) {
  stopifnot(all(c("gene", "avg_log2FC", "p_val_adj") %in% colnames(de_df)))

  # Handle p_val_adj == 0 (or NA) BEFORE computing the rank
  de_df$gene      <- as.character(de_df$gene)
  de_df$p_val_adj <- ifelse(is.na(de_df$p_val_adj) | de_df$p_val_adj <= 0,
                            .Machine$double.xmin, de_df$p_val_adj)
  de_df$rank_stat <- de_df$avg_log2FC * -log10(de_df$p_val_adj)

  # Drop degenerate rows and duplicates, sort descending
  de_df <- de_df[is.finite(de_df$rank_stat) & de_df$rank_stat != 0, ]
  de_df <- de_df[!duplicated(de_df$gene), ]
  de_df <- de_df[order(de_df$rank_stat, decreasing = TRUE), ]

  geneList_symbol <- sort(setNames(de_df$rank_stat, de_df$gene), decreasing = TRUE)

  # SYMBOL -> ENTREZID (gseGO needs Entrez)
  map_tbl <- suppressMessages(
    bitr(names(geneList_symbol), fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db)
  )
  if (is.null(map_tbl) || nrow(map_tbl) == 0) {
    warning("No SYMBOL -> ENTREZ mapping; GSEA returns empty.")
    return(list(ranked_vector = geneList_symbol, result = NULL))
  }
  map_tbl <- dplyr::distinct(tibble::as_tibble(map_tbl), SYMBOL, .keep_all = TRUE)

  geneList_entrez <- geneList_symbol[names(geneList_symbol) %in% map_tbl$SYMBOL]
  names(geneList_entrez) <- map_tbl$ENTREZID[match(names(geneList_entrez), map_tbl$SYMBOL)]
  geneList_entrez <- geneList_entrez[!duplicated(names(geneList_entrez))]
  keep <- !is.na(names(geneList_entrez)) & nzchar(names(geneList_entrez)) &
          is.finite(geneList_entrez)
  geneList_entrez <- sort(geneList_entrez[keep], decreasing = TRUE)

  res <- gseGO(
    geneList      = geneList_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    minGSSize     = minGSSize,
    maxGSSize     = maxGSSize,
    eps           = eps,
    pvalueCutoff  = 1,           # return everything; filter downstream
    pAdjustMethod = pAdjustMethod,
    seed          = seed,
    by            = "fgsea",
    verbose       = verbose
  )

  # Convert Entrez in core_enrichment back to SYMBOL for readability
  if (!is.null(res) && !is.null(res@result) &&
      "core_enrichment" %in% colnames(res@result)) {
    res_df     <- as.data.frame(res@result)
    token_list <- strsplit(res_df$core_enrichment, "/", fixed = TRUE)
    all_tokens <- unique(unlist(token_list))
    all_tokens <- all_tokens[!is.na(all_tokens) & nzchar(all_tokens)]
    if (length(all_tokens) && all(grepl("^[0-9]+$", all_tokens))) {
      id2sym  <- suppressWarnings(
        bitr(all_tokens, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
      )
      map_vec <- setNames(id2sym$SYMBOL, id2sym$ENTREZID)
      res_df$core_enrichment <- sapply(token_list, function(ids) {
        ids  <- ids[!is.na(ids) & nzchar(ids)]
        syms <- map_vec[ids]
        paste(syms[!is.na(syms)], collapse = "/")
      })
      res@result <- res_df
    }
  }
  list(ranked_vector = geneList_entrez, result = res)
}

# ---- make_volcano(): log2FC vs -log10(padj), top genes labelled -------------
make_volcano <- function(de_df, label, outdir,
                         lfc_thr = 1.5, padj_thr = 0.05, n_labels = 25,
                         width = 8, height = 6) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  if (!"gene" %in% colnames(de_df)) de_df <- tibble::rownames_to_column(de_df, "gene")

  df <- de_df
  df$p_val_adj   <- ifelse(is.na(df$p_val_adj) | df$p_val_adj <= 0,
                            .Machine$double.xmin, df$p_val_adj)
  df$nlog10_padj <- -log10(df$p_val_adj)
  df$color_group <- with(df, ifelse(
    avg_log2FC >  lfc_thr & p_val_adj < padj_thr, "Up",
    ifelse(avg_log2FC < -lfc_thr & p_val_adj < padj_thr, "Down", "Other")))

  top_up   <- head(df[df$color_group == "Up", ][
                     order(df[df$color_group == "Up", ]$p_val_adj), ], n_labels)
  top_down <- head(df[df$color_group == "Down", ][
                     order(df[df$color_group == "Down", ]$p_val_adj), ], n_labels)

  p <- ggplot(df, aes(avg_log2FC, nlog10_padj)) +
    geom_point(aes(color = color_group), alpha = 0.6, size = 0.8) +
    geom_vline(xintercept = c(-lfc_thr, 0, lfc_thr),
               linetype = c("dashed", "solid", "dashed"), color = "grey50") +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", color = "grey50") +
    ggrepel::geom_text_repel(data = top_up,   aes(label = gene), size = 3,
                             max.overlaps = Inf, nudge_x = 0.25, hjust = 0) +
    ggrepel::geom_text_repel(data = top_down, aes(label = gene), size = 3,
                             max.overlaps = Inf, nudge_x = -0.25, hjust = 1) +
    scale_color_manual(values = c(Up = "red", Down = "royalblue", Other = "grey70")) +
    labs(title = paste0("Volcano - ", label),
         x = expression(log[2]~"Fold Change"),
         y = expression(-log[10]("adj. p-value"))) +
    theme_minimal(base_size = 11) + theme(legend.position = "none")

  save_plot(p, paste0("Volcano_", clean_name(label)), outdir,
            width = width, height = height, format = "png")
  invisible(p)
}

# ---- plot_gsea_two_panel(): Activated vs Suppressed NES dot plot ------------
.prep_gsea_df <- function(df) {
  req <- c("Description", "NES", "p.adjust", "core_enrichment")
  if (!all(req %in% colnames(df))) stop("GSEA table missing columns.")
  df$direction       <- ifelse(df$NES >= 0, "Activated", "Suppressed")
  df$NES_abs         <- abs(df$NES)
  df$core_gene_count <- ifelse(is.na(df$core_enrichment) | df$core_enrichment == "",
                               0L, lengths(strsplit(df$core_enrichment, "/")))
  df
}

plot_gsea_two_panel <- function(gsea_df, label = "Group", top_n = 12,
                                outdir = ".", width = 10, max_label_chars = 45) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  if (nrow(gsea_df) == 0) { message("No GSEA rows for ", label); return(invisible(NULL)) }
  df <- .prep_gsea_df(gsea_df)

  df_top <- df %>%
    dplyr::group_by(direction) %>%
    dplyr::arrange(p.adjust, dplyr::desc(NES_abs), .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Description = stringr::str_trunc(pretty_label(Description),
                                                   max_label_chars, ellipsis = "...")) %>%
    dplyr::group_by(direction) %>%
    dplyr::mutate(Description_f = forcats::fct_reorder(Description, NES_abs)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(facet_lab = ifelse(direction == "Activated",
                                     paste0("Activated in ", label),
                                     paste0("Suppressed in ", label)))
  if (nrow(df_top) == 0) { message("No enriched terms for ", label); return(invisible(NULL)) }

  n_rows <- max(table(df_top$direction), 1L)
  height <- max(3.5, n_rows * 0.32 + 1.5)

  p <- ggplot(df_top, aes(NES_abs, Description_f)) +
    geom_point(aes(size = core_gene_count, color = p.adjust)) +
    facet_grid(. ~ facet_lab, scales = "free_y") +
    scale_size_continuous(range = c(2, 9), name = "core genes") +
    scale_color_gradient(low = "#e34a33", high = "#3b8bc2",
                         trans = "reverse", name = "p.adjust") +
    labs(title = paste0("GSEA GO BP - ", label),
         x = "Normalized Enrichment Score (|NES|)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey90", color = "grey40"),
          strip.text = element_text(face = "bold"),
          panel.border = element_rect(color = "grey80", fill = NA))

  save_plot(p, paste0("GSEA_twoPanel_", clean_name(label)), outdir,
            width = width, height = height, format = "png")
  invisible(p)
}


# =============================================================================
# WORKFLOW: GATA1s vs wtGATA1, WITHIN each SingleR cell type
# =============================================================================

# ---- 1. Load annotated object -----------------------------------------------
message("\n=== Loading annotated GATA1 object ===")
sobj <- readRDS(io$sobj_rds)
sobj <- JoinLayers(sobj)
DefaultAssay(sobj) <- io$assay
message("  Loaded ", ncol(sobj), " cells")

# Sanity: the columns we key on must exist
stopifnot(io$ident_col %in% colnames(sobj@meta.data),
          io$contrast   %in% colnames(sobj@meta.data))

# Drop the catch-all "other" bucket; test real cell types only
cell_types <- setdiff(unique(as.character(sobj@meta.data[[io$ident_col]])),
                      c("other", NA))
message("  Cell types to test: ", paste(cell_types, collapse = ", "))

# ---- 2. Loop cell types: DE (full ranking) -> GSEA -> plots -----------------
gsea_books <- list()   # collect GSEA tables per cell type
de_books   <- list()   # collect DE tables per cell type

for (ct in cell_types) {
  message("\n--- ", ct, " ---")

  sub <- subset(sobj, cells = colnames(sobj)[sobj@meta.data[[io$ident_col]] == ct])
  grp <- sub@meta.data[[io$contrast]]
  n1  <- sum(grp == io$test_level, na.rm = TRUE)
  n2  <- sum(grp == io$ref_level,  na.rm = TRUE)
  if (n1 < io$min_cells || n2 < io$min_cells) {
    message("  Skipping (", n1, " vs ", n2, " cells - below min_cells).")
    next
  }

  Idents(sub) <- io$contrast
  de <- tryCatch(
    FindMarkers(sub,
                ident.1        = io$test_level,
                ident.2        = io$ref_level,
                test.use       = io$de_test,
                min.pct        = io$min_pct,        # 0 = keep everything for GSEA
                logfc.threshold= io$logfc_thresh,   # 0 = keep everything
                only.pos       = FALSE),
    error = function(e) { message("  DE failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(de) || nrow(de) == 0) next

  de <- tibble::rownames_to_column(de, "gene")
  de_books[[ct]] <- de
  write.csv(de, file.path(outdir_gsea, paste0("DE_", clean_name(ct), ".csv")),
            row.names = FALSE)

  # Volcano from the (unfiltered) DE table
  make_volcano(de, label = paste0(ct, ": ", io$test_level, " vs ", io$ref_level),
               outdir = outdir_gsea,
               lfc_thr = io$volcano_lfc, padj_thr = io$volcano_padj)

  # GSEA
  gs <- run_gsea_go_bp(de, minGSSize = io$gsea_min, maxGSSize = io$gsea_max,
                       eps = io$gsea_eps)
  gsea_df <- to_df(gs$result)
  if (nrow(gsea_df) == 0) { message("  No GSEA terms."); next }

  gsea_books[[ct]] <- gsea_df
  write.csv(gsea_df, file.path(outdir_gsea, paste0("GSEA_", clean_name(ct), ".csv")),
            row.names = FALSE)

  # Two-panel dot plot of the significant terms (padj < 0.05)
  sig <- dplyr::filter(gsea_df, p.adjust < 0.05)
  plot_gsea_two_panel(sig, label = paste0(ct, " (", io$test_level, " vs ", io$ref_level, ")"),
                      outdir = outdir_gsea)
}

# ---- 3. Combined summary across cell types ----------------------------------
if (length(gsea_books)) {
  combined <- dplyr::bind_rows(lapply(names(gsea_books), function(ct) {
    d <- gsea_books[[ct]]; d$cell_type <- ct; d
  }))
  write.csv(combined, file.path(outdir_gsea, "GSEA_all_celltypes.csv"),
            row.names = FALSE)
  message("\nWrote combined GSEA table: ", nrow(combined), " rows across ",
          length(gsea_books), " cell types.")
} else {
  message("\nNo GSEA results produced - check cell counts / mapping.")
}

message("\nOutputs in: ", outdir_gsea)
cat("\nFinished GSEA GO BP workflow.\n")
sessionInfo()
