# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

library('DESeq2')
library('ggplot2')
library('readxl')


# ARGUMENTS
args = strsplit(commandArgs(trailingOnly=TRUE), ' ')[[1]]

root = args[1]
selection_file = args[2]
data_dir = args[3]
hash_value = args[4]
output_dir = args[5]
target_ids = as.numeric(strsplit(args[6], ',')[[1]])
control_ids = as.numeric(strsplit(args[7], ',')[[1]])


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
    s1 = basename(dirname(files[i-1]))
    s2 = basename(dirname(files[i]))
    suffixes = c(paste0('_', s1), paste0('_', s2))
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
    file = paste0(hash_value, '.txt')
    path = file.path(selection_dir, file)
    evaluation_files = c(evaluation_files, path)
  }
}

num_targets = length(target_ids)
num_controls = length(control_ids)

cts = create_table(evaluation_files)
threshold_counts = 10
group_size = floor((num_targets + num_controls) / 2)
filter = rowSums(cts >= threshold_counts) >= group_size
cts = cts[filter, ]

coldata = data.frame(
  condition = c(rep('protein', num_targets), rep('control', num_controls)),
  row.names = colnames(cts)
)
coldata$condition = factor(coldata$condition)

dds = DESeqDataSetFromMatrix(cts, coldata, ~ condition)
dds$condition = relevel(dds$condition, ref = 'control')
dds = DESeq(dds)

resLFC = lfcShrink(dds, coef='condition_protein_vs_control', type='apeglm')
resLFC$targetMean = rowMeans(cts[, 1:num_targets])
resLFC$controlMean = rowMeans(cts[, (num_targets+1):ncol(cts)])
resLFC$enrichFactor = resLFC$targetMean / resLFC$controlMean

resOrdered = resLFC[order(resLFC$padj), ]
resOrdered$negLogp = -log10(resOrdered$padj)
resOrdered$log2FoldChange = resOrdered$log2FoldChange
mask = resOrdered$padj < 0.05
resSign = resOrdered[!is.na(mask) & mask, ]

p = ggplot(resSign, aes(x=log2FoldChange, y=negLogp)) +
  geom_point(alpha=0.5) +
  labs(x='log2 fold change', y='-log10 p-value') +
  xlim(0, ceiling(max(resSign$log2FoldChange)))

cols = c(
  'targetMean',
  'controlMean',
  'enrichFactor',
  'log2FoldChange',
  'padj',
  'negLogp'
)

write.csv(resSign[, cols], file=file.path(output_dir, 'stats.csv'), row.names=TRUE)
ggsave(file.path(output_dir, 'volcano.png'), plot=p, width=15, height=10, dpi=300)

