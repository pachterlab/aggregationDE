library(dplyr)
library(DESeq2)
library(sleuth)
library(tximport)
library(aggregation)

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

files <- sapply(1:6, function(i) file.path(experiment2, exp_string, i, 'kallisto', 'abundance.tsv'))
files <- as.vector(files)
print(files)

design <- c(rep(0,3), rep(1,3))
design <- factor(design)

print('importing tximp on transcript counts')
#tximport Deseq on transcripts
tximp_tx <- tximport(files, type='kallisto', txOut=TRUE)
tx_counts <- round(tximp_tx$counts)
tx_counts <- tx_counts[apply(tx_counts, 1, filter_func), ]
rownames(tx_counts) <- as.character(rownames(tx_counts))
deseq <- DESeq_wt(tx_counts, design)
deseq_results <- as.data.frame(results(deseq))
saveRDS(deseq_results, file.path(experiment, 'wt', exp_string, 'deseq_tx_table.rds'))

deseq_results$target_id <- rownames(deseq_results)
deseq_results <- merge(deseq_results, transcripts, by='target_id', all.x=TRUE)
deseq_results <- mutate(deseq_results, genes = gene_id)
deseq_results <- select(deseq_results, genes, pvalue, baseMean)
deseq_summary <- deseq_results %>% group_by(genes) %>%
	summarise(n = n(),
		tx_deseq = lancaster(pvalue, baseMean),
		tx_deseq_f = fisher(pvalue),
		tx_deseq_min = sidak(pvalue))

files <- sapply(1:6, function(i) file.path(experiment2, exp_string, i, 'kallisto', 'abundance.h5'))
files <- as.vector(files)
print(files)

#sleuth pipeline on transcript counts
sample <- c('s1','s2','s3','s4','s5','s6')
s2c <- data.frame(sample, condition = design, path = files, stringsAsFactors = FALSE)
so <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE) 
sleuth_table <- sleuth_wt_pipeline(so)
if(nrow(sleuth_table)!= nrow(transcripts))
{
	print('error, number of rows in sleuth table and index does not match')
}
sleuth_table <- merge(sleuth_table, transcripts, by='target_id', all.x=TRUE)
saveRDS(sleuth_table, file.path(experiment, 'wt', exp_string, 'sleuth_tx_table.rds'))


sleuth_summary <- sleuth_table %>% group_by(gene_id) %>%
	summarise(tx_sleuth = lancaster(pval,mean_obs),
		tx_sleuth_f = fisher(pval),
		tx_sleuth_min = sidak(pval))
sleuth_summary <- mutate(sleuth_summary, genes = gene_id)
sleuth_summary <- select(sleuth_summary, genes, tx_sleuth, tx_sleuth_f, tx_sleuth_min)

summary_table <- merge(deseq_summary, sleuth_summary, by='genes', all.x=TRUE)

print('finished transcript analysis')
print(dim(summary_table))
saveRDS(summary_table, file.path(experiment, 'wt', exp_string, 'tx_summary_table.rds'))

