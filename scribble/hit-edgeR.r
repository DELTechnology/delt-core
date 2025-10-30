# Auto-generated analysis script
suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
})

args <- list(
  data_path    = "/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/analysis/test/data.csv",
  samples_path = "/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/analysis/test/samples.csv",
  save_dir     = "/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/analysis/test/edgeR",
  cpm          = FALSE
)

# ---- HELPER ----
get_hits.edgeR <- function(data, log = FALSE) {
  
  data.wide <- data %>%
    select(-group) %>%
    pivot_wider(names_from = name, values_from = count, values_fill = 0)
  
  data.row <- data.wide %>%
    select(code_1, code_2)
  
  data.col <- data.frame(
    name = data.wide %>% select(-code_1, -code_2) %>% colnames(),
    stringsAsFactors = FALSE
  )
  groups <- factor(sapply(data.col$name, get_group_from_name))
  groups <- relevel(groups, "naive")
  data.col$group <- groups
  
  data.counts <- as.matrix(data.wide %>% select(-code_1, -code_2))
  rownames(data.counts) <- seq.int(nrow(data.wide))
  
  y <- DGEList(
    counts  = data.counts,
    samples = data.frame(
      name  = data.col$name,
      group = data.col$group,
      row.names = data.col$name
    )
  )
  
  y <- calcNormFactors(y, method = "TMM")
  
  design <- model.matrix(~ 0 + y$samples$group)
  colnames(design) <- levels(y$samples$group)
  
  y   <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  
  cm <- makeContrasts(
    enrichment = protein - no_protein,
    sticky     = no_protein - naive,
    levels = design
  )
  lrt.enrichment <- glmLRT(fit, contrast = cm[, 1])
  lrt.sticky     <- glmLRT(fit, contrast = cm[, 2])
  
  stats.enrichment <- bind_cols(data.row, lrt.enrichment$table) %>%
    mutate(FDR = p.adjust(PValue, method = "BH"))
  stats.sticky <- bind_cols(data.row, lrt.sticky$table) %>%
    mutate(FDR = p.adjust(PValue, method = "BH"))
  
  hits.enrichment <- stats.enrichment %>%
    filter(FDR < 0.05, logFC > 0) %>%
    arrange(desc(logFC), FDR)
  hits.sticky <- stats.sticky %>%
    filter(FDR < 0.05, logFC > 0) %>%
    arrange(desc(logFC), FDR)
  
  counts <- bind_cols(data.row, cpm(y))
  
  list(
    stats  = list(enrichment = stats.enrichment,
                  sticky     = stats.sticky),
    hits   = list(enrichment = hits.enrichment,
                  sticky     = hits.sticky),
    counts = counts
  )
}

# ---- Load data ----
data    <- read_csv(args$data_path, show_col_types = FALSE)
samples <- read_csv(args$samples_path, show_col_types = FALSE)

grp_by_name <- setNames(samples$group, samples$name)
get_group_from_name <- function(name) grp_by_name[[name]]

data <- data %>%
  inner_join(samples, by = "name")

# ---- ANALYSIS ----
result.edgeR <- get_hits.edgeR(data = data)

dir.create(args$save_dir, showWarnings = FALSE, recursive = TRUE)

# save stats
for (i in seq_along(result.edgeR$stats)) {
  name  <- names(result.edgeR$stats)[i]
  stats <- result.edgeR$stats[[i]]
  save.path <- file.path(args$save_dir, paste0(name, "_stats.csv"))
  write_csv(stats, file = save.path)
}

# save hits
for (i in seq_along(result.edgeR$hits)) {
  name <- names(result.edgeR$hits)[i]
  hits <- result.edgeR$hits[[i]]
  save.path <- file.path(args$save_dir, paste0(name, "_hits.csv"))
  write_csv(hits, file = save.path)
}

# derive selections from result table columns and export per-selection counts
selections <- setdiff(colnames(result.edgeR$counts), c("code_1", "code_2"))

# save normalized counts
for (selection in selections) {
  fname <- paste0(selection, ".csv")
  save.path <- file.path(args$save_dir, fname)
  result.edgeR$counts %>%
    select(code_1, code_2, all_of(selection)) %>%
    mutate(count = .data[[selection]]) %>%
    select(-all_of(selection)) %>%
    write_csv(file = save.path)
}