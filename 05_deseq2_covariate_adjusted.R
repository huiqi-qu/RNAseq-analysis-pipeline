#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

counts_path <- "raw_counts.csv"
metadata_path <- "metadata.csv"

count_data <- fread(counts_path)
count_mat <- as.matrix(
  data.frame(count_data, row.names = 1, check.names = FALSE)
)
mode(count_mat) <- "numeric"
count_mat <- round(count_mat)

keep <- rowSums(count_mat) > 1
count_mat <- count_mat[keep, , drop = FALSE]

meta <- read.csv(metadata_path, stringsAsFactors = FALSE)

required <- c("Sample", "batch", "AGE", "Sex", "Race", "BigGroups")
missing <- setdiff(required, names(meta))
if (length(missing) > 0) {
  stop("Missing metadata columns: ", paste(missing, collapse = ", "))
}

rownames(meta) <- meta$Sample

meta$AGE <- suppressWarnings(as.numeric(meta$AGE))

meta$Sex <- tolower(trimws(meta$Sex))
meta$Sex[meta$Sex %in% c("", "na", "n/a")] <- NA
meta$Sex <- factor(meta$Sex, levels = c("male", "female"))

meta$Race <- trimws(meta$Race)
meta$Race[meta$Race == "w"] <- "W"
meta$Race[meta$Race %in% c("", "NA", "na", "N/A")] <- NA
meta$Race[!(meta$Race %in% c("B", "W"))] <- NA
meta$Race <- factor(meta$Race, levels = c("W", "B"))

meta$batch <- factor(meta$batch)
meta$BigGroups <- factor(
  meta$BigGroups,
  levels = c("Control", "Ischemic", "Nonischemic")
)

common <- intersect(colnames(count_mat), rownames(meta))
count_mat <- count_mat[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]

design_cols <- c("batch", "AGE", "Sex", "Race", "BigGroups")
complete <- complete.cases(meta[, design_cols])

if (sum(!complete) > 0) {
  message("Excluding ", sum(!complete), " samples with missing design covariates.")
}

meta <- meta[complete, , drop = FALSE]
count_mat <- count_mat[, rownames(meta), drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta,
  design = ~ batch + AGE + Sex + Race + BigGroups
)

dds <- DESeq(dds)

res_is_ctrl <- results(
  dds,
  contrast = c("BigGroups", "Ischemic", "Control")
)
res_ni_ctrl <- results(
  dds,
  contrast = c("BigGroups", "Nonischemic", "Control")
)
res_ni_is <- results(
  dds,
  contrast = c("BigGroups", "Nonischemic", "Ischemic")
)

norm_counts <- counts(dds, normalized = TRUE)
mean_expression <- rowMeans(norm_counts, na.rm = TRUE)

res_is_ctrl$mean_expression <- mean_expression
res_ni_ctrl$mean_expression <- mean_expression
res_ni_is$mean_expression <- mean_expression

res_is_ctrl <- res_is_ctrl[order(res_is_ctrl$padj), ]
res_ni_ctrl <- res_ni_ctrl[order(res_ni_ctrl$padj), ]
res_ni_is <- res_ni_is[order(res_ni_is$padj), ]

write.csv(
  as.data.frame(res_is_ctrl),
  "DE_Ischemic_vs_Control_covariate_adjusted.csv"
)
write.csv(
  as.data.frame(res_ni_ctrl),
  "DE_Nonischemic_vs_Control_covariate_adjusted.csv"
)
write.csv(
  as.data.frame(res_ni_is),
  "DE_Nonischemic_vs_Ischemic_covariate_adjusted.csv"
)

vst_obj <- vst(dds, blind = TRUE)
write.csv(
  assay(vst_obj),
  "VST_counts_for_QC.csv",
  quote = FALSE
)

message("Covariate-adjusted DESeq2 analysis complete.")
