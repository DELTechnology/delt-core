library(tidyverse)
library(DESeq2)
library(edgeR)
library(data.table)

# %%
data.dir <- file.path("/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/selections")
save.dir = file.path("/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/analysis/")
dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)

selections = c("AG24_1", "AG24_2", "AG24_3")
counts = read.delim(file.path(data.dir, "AG24_1", "counts.txt"), sep = "\t")

# %%
read_counts <- function(name) {
  df <- read.delim(file.path(sel_dir, name, "counts.txt"), sep = "\t")
  df$name <- name
  df
}

get_cpm = function(counts) {
  counts$count = counts$count / sum(counts$count) * 1e6
  counts
}

get_deseq2 = function(selections) {
  # NOTE: collapseReplicates can be used if you have technical replicates
  #   how should we treat our replicas?
  # DOCS: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat
  
  data = bind_rows(
    lapply(selections, read_counts),
  )
  
  data = data %>% 
    pivot_wider(names_from = name, values_from = count, values_fill = 0)
  
  data.row = data %>% 
    select(code_1, code_2)
  
  data.col = data.frame(
    row.names = data %>% select(-code_1, -code_2) %>% colnames() %>% factor(),
    group = factor(rep("A", length(selections)))
    )
  
  cm = as.matrix(data %>% select(-code_1, -code_2))
    
  dds <- DESeqDataSetFromMatrix(
    countData = cm,
    colData = data.col,
    design = ~ 1)
  
  dds <- DESeq2::estimateSizeFactors(dds)
  data = counts(dds, normalized = TRUE)
  colnames(data) = colData(dds) %>% rownames()
  bind_cols(data.row, as_tibble(data, rownames = NA))
}

get_tmm = function(selections, log=FALSE) {
  
  data = bind_rows(
    lapply(selections, read_counts),
  )
  
  data = data %>% 
    pivot_wider(names_from = name, values_from = count, values_fill = 0)
  
  data.row = data %>% 
    select(code_1, code_2)
  
  data.col = data.frame(
    row.names = data %>% select(-code_1, -code_2) %>% colnames() %>% factor(),
    group = factor(rep("A", length(selections)))
  )
  
  cm = as.matrix(data %>% select(-code_1, -code_2))
  
  y <- DGEList(counts = cm)
  y <- calcNormFactors(y, method = "TMM")
  data = edgeR::cpm(y, normalized.lib.sizes = TRUE, log = log, prior.count = prior.count)
  
  colnames(data) = colData(dds) %>% rownames()
  bind_cols(data.row, as_tibble(data, rownames = NA))
}

# %% CPM
for(selection in selections){
  counts = read.delim(file.path(data.dir, selection, "counts.txt"), sep = "\t")
  cpm = get_cpm(counts)
  save.path = file.path(data.dir, selection, "cpm.txt")
  write_delim(cpm, file = save.path, delim = '\t')
}

bb.names = c('code_1', 'code_2')
data = get_deseq2(selections)
for(selection in selections){
  
  d = as.data.table(data)
  
  cols = c(bb.names, selection)
  d = d[, ..cols]
  setnames(d, selection, "count")
  
  save.path = file.path(data.dir, selection, "deseq2.txt")
  fwrite(d, file = save.path, sep = "\t")
}

data = get_tmm(selections)
for(selection in selections){
  
  d = as.data.table(data)
  
  cols = c(bb.names, selection)
  d = d[, ..cols]
  setnames(d, selection, "count")
  
  save.path = file.path(data.dir, selection, "tmm.txt")
  fwrite(d, file = save.path, sep = "\t")
}
