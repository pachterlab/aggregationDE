library(dplyr)
library(DESeq2)
library(sleuth)
library(limma)
library(tximport)

source('misc.R')

#experiment <- '~/geuvadis/gfr_3_3_20_42_2'
#experiment2 <- '/home/hjp/sleuth_paper_analysis/geuvadis/sims/gfr_3_3_20_42_2'
#exp_string <- 'exp_1'
#exp_num <- 1
#
args <- commandArgs(TRUE)
experiment <- args[[1]]
experiment2 <- args[[2]]
exp_string <- args[[3]]
exp_num <- as.numeric(args[[4]])
print(experiment)
print(experiment2)
print(exp_string)
print(exp_num)

# read transcript to gene map
# one to one table
transcripts <- read.table('~/indices/Homo_sapiens.GRCh38.rel79.transcripts', sep='\t', stringsAsFactors = FALSE)
colnames(transcripts) <- c('target_id', 'gene_id')

design <- c(rep(0,3), rep(1,3))
design <- factor(design)

print('starting gene pipeline')

#tximport Gene Counts
print('importing kallisto files as gene counts')
files <- sapply(1:6, function(i) file.path(experiment2, exp_string, i, 'kallisto', 'abundance.tsv'))
files <- as.vector(files)
print(files)
tximp <- tximport(files, type='kallisto', tx2gene=transcripts)
gene_counts <- round(tximp$counts)
gene_counts <- gene_counts[apply(gene_counts, 1, filter_func), ]

#run DESeq on Gene Counts
deseq <- DESeq_wt(gene_counts, design)
deseq_results <- as.data.frame(results(deseq))
saveRDS(deseq_results, file.path(experiment, 'wt', exp_string, 'deseq_gene_table.rds'))

deseq_results$genes <- rownames(deseq_results)
deseq_results <- mutate(deseq_results, gene_deseq = pvalue)
deseq_results <- select(deseq_results, genes, gene_deseq)


print('starting sleuth gene pipeline')
files <- sapply(1:6, function(i) file.path(experiment2, exp_string, i, 'kallisto', 'abundance.h5'))
files <- as.vector(files)
print(files)
sample <- c('s1','s2','s3','s4','s5','s6')
s2c <- data.frame(sample, condition = design, path = files, stringsAsFactors = FALSE)
so <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE, target_mapping = transcripts, aggregation_column = 'gene_id')
sleuth_table <- sleuth_wt_pipeline(so)
saveRDS(sleuth_table, file.path(experiment, 'wt', exp_string, 'sleuth_gene_table.rds'))

sleuth_table <- mutate(sleuth_table, genes = target_id, gene_sleuth = pval)
sleuth_table <- select(sleuth_table, genes, gene_sleuth)
summary_table <- merge(deseq_results, sleuth_table, by='genes', all.x=TRUE)
saveRDS(summary_table, file.path(experiment, 'wt', exp_string, 'gene_summary_table.rds'))

print('finished gene analysis')
print(dim(summary_table))

