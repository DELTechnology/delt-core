library('DESeq2')

setwd('/Users/Gary/Downloads/zivi/data/normalization')
file = 'data/counts.txt'
cts = as.matrix(read.csv(file, sep='\t', header=TRUE, row.names='Code'))

coldata = data.frame(
  condition = c('protein', 'protein', 'protein', 'control', 'control', 'control'),
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
resultsNames(dds)
resLFC = lfcShrink(dds, coef='condition_protein_vs_control', type='apeglm')

resOrdered = resLFC[order(resLFC$pvalue), ]
sum(resOrdered$padj < 0.1, na.rm=TRUE)
head(as.data.frame(resOrdered))
