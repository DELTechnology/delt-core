# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

library('DESeq2')
library('readxl')
library('tidyverse')


# ARGUMENTS
args = strsplit(commandArgs(trailingOnly=TRUE), ' ')[[1]]

root = args[1]
selection_file = args[2]
data_dir = args[3]
target_ids = as.numeric(strsplit(args[4], ',')[[1]])
control_ids = as.numeric(strsplit(args[5], ',')[[1]])


# FUNCTIONS
read_selection = function(file) {
  selection = read.table(file, header=TRUE, sep='\t')
  num_cols = dim(selection)[2]
  selection$Code = selection$Code1
  for (col in 3:(num_cols)) {
    selection$Code = paste(selection$Code, selection[, col], sep='-')
  }
  return(selection[, c('Count', 'Code')])
}

create_table = function(files) {
  selections = data.frame(Code=character())
  for (i in seq_along(files)) {
    selection = read_selection(files[i])
    suffixes = c(paste('_', files[i-1], sep=''), paste('_', files[i], sep=''))
    selections = merge(selections, selection, by='Code', all.x=TRUE, all.y=TRUE, suffixes=suffixes)
  }
  selections[is.na(selections)] = 0
  rownames(selections) = selections$Code
  return(selections[, -1])
}


# MAIN
setwd(root)
selections = read_excel(selection_file)
selection_ids = c(target_ids, control_ids)
evaluation_files = c()

for (condition in selection_ids) {
  for (id in condition) {
    selection_dir = file.path(data_dir, paste0('selection-', id))
    files = list.files(selection_dir)
    file = files[grepl('\\.txt$', files)]
    path = file.path(selection_dir, file)
    evaluation_files = c(evaluation_files, path)
  }
}

num_targets = length(target_ids)
num_controls = length(control_ids)

cts = create_table(evaluation_files)
coldata = data.frame(
  condition = c(rep('protein', num_targets), rep('control', num_controls)),
  row.names = colnames(cts)
)
coldata$condition = factor(coldata$condition)

dds = DESeqDataSetFromMatrix(cts, coldata, ~ condition)

threshold_counts = 10
group_size = floor((num_targets + num_controls) / 2)
filter = rowSums(counts(dds) >= threshold_counts) >= group_size
dds = dds[filter, ]
dds$condition = relevel(dds$condition, ref = 'control')

dds = DESeq(dds)

res = results(dds, name='condition_protein_vs_control')
resLFC = lfcShrink(dds, coef='condition_protein_vs_control', type='apeglm')
resOrdered = resLFC[order(resLFC$padj), ]
resOrdered$neglogp = -log10(resOrdered$padj)
mask = resOrdered$padj < 0.05
resSign = resOrdered[!is.na(mask) & mask, ]

p = ggplot(resOrdered, aes(x=log2FoldChange, y=neglogp)) +
  geom_point(alpha=0.5) +
  labs(x='log2 fold change', y='-log10 p-value')

write.csv(resSign, file='stats.csv', row.names=TRUE)
ggsave('volcano.png', plot=p, width=6, height=4, dpi=300)

