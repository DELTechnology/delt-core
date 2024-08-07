library('DESeq2')
library('tidyverse')

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
path = '/Volumes/data_shared/data/shared/AM_Andy_DataAnalysisQuatsch/Data Campaign AG_21'
setwd(path)

files = c(
  'selection_12_.txt', # protein
  'selection_13_.txt', # protein
  'selection_14_.txt', # protein
  'selection_8_.txt',  # control
  'selection_19_.txt', # control
  'selection_20_.txt'  # control
)

cts = create_table(files)

coldata = data.frame(
  condition = c(rep('protein', 3), rep('control', 3)),
  row.names = colnames(cts)
)
coldata$condition = factor(coldata$condition)

dds = DESeqDataSetFromMatrix(cts, coldata, ~ condition)

group_size = 3
filter = rowSums(counts(dds) >= 10) >= group_size
dds = dds[filter, ]
dds$condition = relevel(dds$condition, ref = 'control')

dds = DESeq(dds)

res = results(dds, name='condition_protein_vs_control')
# plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
resLFC = lfcShrink(dds, coef='condition_protein_vs_control', type='apeglm')

resOrdered = resLFC[order(resLFC$padj), ]
resOrdered$neglogp = -log10(resOrdered$padj)
mask = resOrdered$padj < 0.05
resSign = resOrdered[!is.na(mask) & mask, ]

head(as.data.frame(resSign))
# write.csv(resSign, file = "stats.csv", row.names = TRUE)

ggplot(resOrdered, aes(x=log2FoldChange, y=neglogp)) +
  geom_point(alpha=0.5) +
  labs(
    x = 'log2 fold change',
    y = '-log10 p-value'
  ) +
  ylim(0, 5) +
  theme_minimal()



# this doesn't work

files = c(
  'selection_12_.txt', # protein
  'selection_8_.txt'   # control
)

cts = create_table(files)

coldata = data.frame(
  condition = c(rep('protein', 1), rep('control', 1)),
  row.names = colnames(cts)
)
coldata$condition = factor(coldata$condition)

dds = DESeqDataSetFromMatrix(cts, coldata, ~ condition)

group_size = 1
filter = rowSums(counts(dds) >= 10) >= group_size
dds = dds[filter, ]
dds$condition = relevel(dds$condition, ref = 'control')

dds = DESeq(dds)

res = results(dds, name='condition_protein_vs_control')
# plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
resLFC = lfcShrink(dds, coef='condition_protein_vs_control', type='apeglm')

resOrdered = resLFC[order(resLFC$padj), ]
resOrdered$neglogp = -log10(resOrdered$padj)
mask = resOrdered$padj < 0.05
resSign = resOrdered[!is.na(mask) & mask, ]

head(as.data.frame(resSign))
# write.csv(resSign, file = "stats.csv", row.names = TRUE)

ggplot(resOrdered, aes(x=log2FoldChange, y=neglogp)) +
  geom_point(alpha=0.5) +
  labs(
    x = 'log2 fold change',
    y = '-log10 p-value'
  ) +
  ylim(0, 5) +
  theme_minimal()

