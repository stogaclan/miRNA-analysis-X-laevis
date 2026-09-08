# ---------------------------------------------------------------------------
# Differential expression of X. laevis miRNAs — spinal cord regeneration
#
# Input:  mirdeep2/output/mirna_count_matrix.tsv  (from merge_quantifier_outputs.py)
#         metadata.csv
# Output: deseq2_results/*.csv, *.pdf
#
# Two tissues (NPSC, fluid) analysed separately, three contrasts each.
# EXPLORATORY: n = 2 per group, thresholded on raw p. See section 12 of index.md.
#
# Work through index.md the first time. Use this to re-run.
# ---------------------------------------------------------------------------

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# --- EDIT THESE ------------------------------------------------------------
count_file <- "mirdeep2/output/mirna_count_matrix.tsv"
min_count  <- 10     # a miRNA needs this many counts...
min_samps  <- 2      # ...in at least this many samples (= one full group)
p_cut      <- 0.05   # RAW p, not padj — justified in index.md section 12
lfc_cut    <- 1      # log2(2) = 2-fold
# ---------------------------------------------------------------------------

dir.create("deseq2_results", recursive = TRUE, showWarnings = FALSE)

# --- Load ------------------------------------------------------------------
raw  <- read.delim(count_file, check.names = FALSE)

# strip.white + trimws: stray spaces after the commas in metadata.csv would
# otherwise make "npsc" and " npsc" different values, and nothing would match.
meta <- read.csv("metadata.csv", row.names = 1, strip.white = TRUE)
meta[] <- lapply(meta, function(x) if (is.character(x)) trimws(x) else x)
rownames(meta) <- trimws(rownames(meta))

meta$condition <- factor(meta$condition, levels = c("uninjured", "dpt2", "dpt6"))
meta$tissue    <- factor(meta$tissue)

# merge_quantifier_outputs.py has already summed each miRNA's precursor rows
# and rounded the fractional -W counts, so these are ready to use.
missing <- setdiff(rownames(meta), colnames(raw))
if (length(missing))
  stop("Samples in metadata.csv but not in the count matrix: ",
       paste(missing, collapse = ", "),
       "\n  Matrix columns are: ", paste(setdiff(colnames(raw), "miRNA"), collapse = ", "))

counts_all <- as.matrix(raw[, rownames(meta)])
rownames(counts_all) <- raw$miRNA
message(nrow(counts_all), " miRNAs x ", ncol(counts_all), " samples")

# --- Fit one model per tissue ----------------------------------------------
run_tissue <- function(tis) {
  m <- meta[meta$tissue == tis, , drop = FALSE]

  if (nrow(m) < 2)
    stop("Found ", nrow(m), " samples with tissue == '", tis,
         "'. metadata.csv contains: ",
         paste(sort(unique(as.character(meta$tissue))), collapse = ", "),
         "\n  Check spelling and capitalisation in the tissue column.")

  m$condition <- droplevels(m$condition)
  cts <- counts_all[, rownames(m), drop = FALSE]

  dds  <- DESeqDataSetFromMatrix(cts, m, design = ~ condition)
  keep <- rowSums(counts(dds) >= min_count) >= min_samps
  dds  <- dds[keep, ]
  message(tis, ": ", nrow(dds), " miRNAs kept after filtering")

  DESeq(dds)
}

contrasts <- list(
  dpt2_vs_uninjured = c("condition", "dpt2", "uninjured"),
  dpt6_vs_uninjured = c("condition", "dpt6", "uninjured"),
  dpt6_vs_dpt2      = c("condition", "dpt6", "dpt2")
)

analyse <- function(tis) {

  dds <- run_tissue(tis)

  # vst() needs >= 1000 features by default and errors on a few hundred miRNAs.
  # varianceStabilizingTransformation() is the same transform without that limit.
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

  pdf(sprintf("deseq2_results/%s_qc.pdf", tis), width = 7, height = 6)
  print(plotPCA(vsd, intgroup = "condition") + ggtitle(paste(tis, "- condition")))
  plotDispEsts(dds)
  dev.off()

  for (nm in names(contrasts)) {
    r <- as.data.frame(results(dds, contrast = contrasts[[nm]]))
    r$miRNA <- rownames(r)
    r <- r[order(r$pvalue), c("miRNA", "baseMean", "log2FoldChange",
                              "lfcSE", "pvalue", "padj")]

    # RAW p. padj is kept in the full table — check it before making claims.
    sig <- subset(r, !is.na(pvalue) & pvalue < p_cut & abs(log2FoldChange) > lfc_cut)

    write.csv(r,   sprintf("deseq2_results/%s_%s_all.csv", tis, nm), row.names = FALSE)
    write.csv(sig, sprintf("deseq2_results/%s_%s_sig.csv", tis, nm), row.names = FALSE)

    n_padj <- sum(!is.na(r$padj) & r$padj < 0.05)
    message(tis, " ", nm, ": ", nrow(sig), " candidates (up ",
            sum(sig$log2FoldChange > 0), ", down ",
            sum(sig$log2FoldChange < 0), ")",
            if (n_padj) paste0("  [", n_padj, " also survive padj < 0.05]") else "")

    # volcano
    r$hit <- !is.na(r$pvalue) & r$pvalue < p_cut & abs(r$log2FoldChange) > lfc_cut
    p <- ggplot(r, aes(log2FoldChange, -log10(pvalue))) +
      geom_point(aes(colour = hit), alpha = 0.7) +
      geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
      geom_hline(yintercept = -log10(p_cut),        linetype = "dashed") +
      geom_text_repel(data = subset(r, hit), aes(label = miRNA),
                      size = 2.5, max.overlaps = 20) +
      scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
      labs(x = "log2 fold change", y = "-log10 raw p",
           title    = paste(tis, nm),
           subtitle = "Exploratory, n = 2 per group; raw p-values") +
      theme_bw() + theme(legend.position = "none")
    ggsave(sprintf("deseq2_results/%s_%s_volcano.pdf", tis, nm), p,
           width = 6, height = 5)

    # heatmap of top candidates
    if (nrow(sig) >= 2) {
      top <- head(sig$miRNA, 30)
      ann <- meta[colnames(vsd), "condition", drop = FALSE]
      pdf(sprintf("deseq2_results/%s_%s_heatmap.pdf", tis, nm), width = 7, height = 8)
      pheatmap(assay(vsd)[top, ], scale = "row", annotation_col = ann,
               main = paste(tis, nm))
      dev.off()
    }
  }
  dds
}

dds_npsc  <- analyse("npsc")
dds_fluid <- analyse("fluid")

message("\nDone. See deseq2_results/")
message("Reminder: these are CANDIDATES, not differentially expressed miRNAs. ",
        "n = 2 per group, thresholded on raw p. Validate before believing.")
