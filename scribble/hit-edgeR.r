# %%
library(tidyverse)
library(edgeR)
library(data.table)

# %%
name = 'test'
data.dir <- file.path("/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/selections")
save.dir = file.path("/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/analysis/", name, "edgeR")
dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)

selections = c("AG24_1", "AG24_2", "AG24_3",
               "AG24_10", "AG24_11", "AG24_12",
               "AG24_19", "AG24_20", "AG24_21"
               )
grps = groups = c(
  rep('no_protein', 3),
  rep('protein', 3),
  rep('naive', 3)
  )

get_group_from_name = function(name){
  idx = which(selections == name)
  grps[idx]
}

# %%
read_counts <- function(name) {
  df <- read.delim(file.path(data.dir, name, "counts.txt"), sep = "\t")
  df$name <- name
  df
}

read_dataset = function(selections, groups){
  data <- bind_rows(
    Map(function(sel, grp) {
      d <- read_counts(sel)
      d$group <- grp
      d
    }, selections, groups)
  )
  
  data$group <- factor(data$group)
  data
}

get_hits.edgeR = function(data, log=FALSE) {
  
  data.wide = data %>% 
    select(-group) %>%
    pivot_wider(names_from = name, values_from = count, values_fill = 0)
  
  data.row = data.wide %>% 
    select(code_1, code_2)
  # data.row$feature_id <- seq_len(nrow(data.row))
  
  data.col = data.frame(
    name = data.wide %>% select(-code_1, -code_2) %>% colnames() %>% factor()
  )
  groups = factor(sapply(data.col$name, get_group_from_name))
  groups = relevel(groups, "naive")
  data.col$group = groups
  
  data.counts = as.matrix(data.wide %>% select(-code_1, -code_2))
  rownames(data.counts) = seq.int(nrow(data.wide))
  
  y <- DGEList(counts = data.counts,
               samples = data.col$name,
               group = data.col$group)
  y <- calcNormFactors(y, method = "TMM")
  
  design <- model.matrix(~ 0 + y$samples$group)
  colnames(design) <- levels(y$samples$group)
  
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  
  cm = makeContrasts(
    enrichment = protein - no_protein,
    sticky = no_protein - naive,
    levels = design
  )
  lrt.enrichment <- glmLRT(fit, contrast = cm[, 1])
  lrt.sticky <- glmLRT(fit, contrast = cm[, 2])
  
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
  
  counts = bind_cols(data.row, edgeR::cpm(y))
  
  result = list(
    stats = list(enrichment=stats.enrichment,
                 sticky=stats.sticky),
    hits = list(enrichment=hits.enrichment,
                sticky=hits.sticky),
    counts = counts
  )
  
}

data = read_dataset(selections, grps)
result.edgeR = get_hits.edgeR(data=data)

for (i in seq_along(result.edgeR$stats)) {
  name  <- names(result.edgeR$stats)[i]
  stats <- result.edgeR$stats[[i]]
  
  save.path <- file.path(save.dir, paste0(name, "_stats.csv"))
  write_csv(stats, file = save.path)
}

for (i in seq_along(result.edgeR$hits)) {
  name  <- names(result.edgeR$hits)[i]
  hits <- result.edgeR$hits[[i]]
  
  save.path <- file.path(save.dir, paste0(name, "_hits.csv"))
  write_csv(hits, file = save.path)
}

for(selection in selections){
  fname = paste0(selection, '.csv')
  save.path = file.path(save.dir, fname)
  result.edgeR$counts %>% 
    select(code_1, code_2, all_of(selection)) %>%
    mutate(count = selection) %>% 
    select(-all_of(selection)) %>%
    write_csv(file=save.path)
}

