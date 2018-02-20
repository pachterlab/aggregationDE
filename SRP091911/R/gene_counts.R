library(sleuth)
library(aggregation)
source('~/TCC/misc.R')

base_path <- '/home/lynnyi/SRP091911'

print('reading ec map')
ec_map_path <- file.path(base_path, 'pseudoquant/matrix.ec')
print(ec_map_path)
ec_map <- read.table(ec_map_path, sep='\t', stringsAsFactors=FALSE)
colnames(ec_map) <- c('tcc_num', 'transcript_num')

print('reading matrix')
in_matrix <- file.path(base_path, 'pseudoquant/matrix.tsv')
print(in_matrix)
df <- scan(in_matrix, skip = 1, what=numeric(), sep='\t')
df <- matrix(df, nrow=21)
df <- df[-1,]
print(dim(df))
print('first 10 rows/samples')
df <- df[c(1:10),]

print('reading transcript_gene map')
transcripts <- read.table('/home/lynnyi/transcriptomes/Rattus_norvegcius.Rnor_6.0.transcripts', stringsAsFactors=FALSE)
colnames(transcripts) <- c('target_id', 'genes')
transcripts$target_id <- gsub('\\..*', '', transcripts$target_id)
transcripts$genes <- gsub('\\..*', '', transcripts$genes)

tcc_table <- readRDS(file.path(base_path, 'wt', 'tcc_table.rds'))
lan_table <- readRDS(file.path(base_path,'wt','tcc_pipeline_results.rds'))

colnames(df) <- tcc_table$genes
counts <- t(rowsum(t(df), colnames(df)))
counts <- counts[,-1]

design <- c(rep(0, 5), rep(1,5))
design <- factor(design)
deseq <- DESeq_wt(t(counts), design)
deseq_results <- results(deseq)

results <- data.frame(gene = rownames(deseq_results), gene_pval = deseq_results$pval, lan_pval = lan_table$pval)

png(file.path(base_path, 'wt', 'plots', 'pval_scatter.png'), height = 5000, width = 5000, res = 1200)
p <- ggplot(data = results, aes(x=-log(lan_pval), y=-log(gene_pval))) + geom_point(alpha = .1) + xlab('Lancaster TCC (-log(p-value))') + ylab('Gene Counts (-log(p-value))')
print(p)
dev.off()


