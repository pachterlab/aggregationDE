library(dplyr)
library(DESeq2)
library(sleuth)
library(limma)
library(tximport)

source('aggregation.R')
source('roc_curve.R')

#experiment <- '~/geuvadis/gfr_3_3_20_42_2'
#experiment2 <- '/home/hjp/sleuth_paper_analysis/geuvadis/sims/gfr_3_3_20_42_2'
#exp_string <- 'exp_1'
#exp_num <- 1

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

# read matrix.ec
# make sure to read as strings
# one to many table
# this one is 0 indexed
ec_map_path <- file.path(experiment, exp_string, 'matrix.ec')
print('reading ec map')
print(ec_map_path)
ec_map <- read.table(ec_map_path, sep='\t', stringsAsFactors=FALSE)
colnames(ec_map) <- c('tcc_num', 'transcript_num')

in_matrix <- file.path(experiment, exp_string, 'matrix.tsv')
print('reading matrix')
print(in_matrix)
df <- scan(in_matrix, skip = 1, what=numeric(), sep='\t')
df <- matrix(df, nrow=7)
df <- df[-1,]
print(dim(df))

# call DESeq on TCCs
design <- c(rep(0,3), rep(1,3))
design <- factor(design)

filter_func <- function(row)
{
	mean(row) > 0
}

# need ncol counts == nrow design == nsamples
#print('starting DESeq on TCC analysis')
DESeqData <- DESeq(t(df), design)
dresults <- as.data.frame(results(DESeqData))
deseq_pvalues <- dresults$pvalue

print('starting limma on TCC analysis')
lresults <- limma(t(df), design)
limma_pvalues <- lresults$p.value[,2]

#sleuth_tcc
print('starting sleuth on TCC analysis')
tcc_files <- sapply(1:6, function(i) file.path(experiment, exp_string, 'sleuth_tcc', i, 'abundance.h5'))
print(tcc_files)
sample <- c('s1','s2','s3','s4','s5','s6')
s2c <- data.frame(sample, condition = design, path = tcc_files, stringsAsFactors = FALSE)
print(s2c)
so <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE, filter_fun = filter_func)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_wt(so, 'condition1', 'full')

sleuth_table <- sleuth_results(so, 'condition1', which_model = 'full')
sleuth_table$target_id <- as.numeric(sleuth_table$target_id)
sleuth_table <- sleuth_table[order(sleuth_table$target_id),] 
sleuth_pvalues <- sleuth_table$pval

#map tccs to genes
tcc2genemap <- make_tcc_to_gene_map2(transcripts, ec_map)
tcc_table <- data.frame(genes = tcc2genemap, dpvalue = deseq_pvalues, lpvalue = limma_pvalues, spvalue = sleuth_pvalues, meancounts = colMeans(df))

summary_table <- tcc_table %>% group_by(genes) %>%
	summarise(n = n(),
		tcc_deseq = lancaster(dpvalue, meancounts),
		tcc_limma = lancaster(lpvalue, meancounts),
		tcc_sleuth = lancaster(spvalue, meancounts))
summary_table <- summary_table[-1,]

print('finished TCC analyses')
print(head(summary_table))
saveRDS(summary_table, file.path(experiment, exp_string, 'summary_table1.rds'))
saveRDS(tcc_table, file.path(experiment, exp_string, 'tcc_table.rds'))


print('starting gene pipeline')

#tximport Gene Counts
print('importing kallisto files as gene counts')
files <- sapply(1:6, function(i) file.path(experiment2, exp_string, i, 'kallisto', 'abundance.h5'))
files <- as.vector(files)
print(files)
tximp <- tximport(files, type='kallisto', tx2gene=transcripts)
gene_counts <- round(tximp$counts)
gene_counts <- gene_counts[apply(gene_counts, 1, filter_func), ]

#run DESeq on Gene Counts
tximp_deseq <- DESeq(gene_counts, design)
tximp_results <- as.data.frame(results(tximp_deseq))
tximp_results$genes <- rownames(tximp_results)
tximp_results <- mutate(tximp_results, gene_deseq = pvalue)
tximp_results <- select(tximp_results, genes, gene_deseq)
#summary_table$deseq is now the pvalues for this pipeline
summary_table <- merge(summary_table, tximp_results, by = 'genes', all.x= TRUE)
#run limma on Gene Counts
tximp_limma <- limma(gene_counts, design)
tximp_limma_results <- data.frame(gene_limma = tximp_limma$p.value[,2], genes = as.character(rownames(gene_counts)))
summary_table <- merge(summary_table, tximp_limma_results, by='genes', all.x=TRUE)


print('starting sleuth gene pipeline')
s2c <- data.frame(sample, condition = design, path = files, stringsAsFactors = FALSE)
so <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE, target_mapping = transcripts, aggregation_column = 'gene_id')
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_wt(so, 'condition1', 'full')
sleuth_table <- sleuth_results(so, 'condition1', which_model = 'full')
sleuth_table <- mutate(sleuth_table, genes = target_id, gene_sleuth = pval)
sleuth_table <- select(sleuth_table, genes, gene_sleuth)
summary_table <- merge(summary_table, sleuth_table, by='genes', all.x=TRUE)

print('finished gene analysis')
print(dim(summary_table))


print('importing tximp on transcript counts')
#tximport Deseq on transcripts
tximp_tx <- tximport(files, type='kallisto', txOut=TRUE)
tx_counts <- round(tximp_tx$counts)
tx_counts <- tx_counts[apply(tx_counts, 1, filter_func), ]
rownames(tx_counts) <- as.character(rownames(tx_counts))
tximp_deseq_tx <- DESeq(tx_counts, design)
tximp_tx_results <- as.data.frame(results(tximp_deseq_tx))
tximp_tx_results$target_id <- rownames(tximp_tx_results)
tximp_tx_results <- merge(tximp_tx_results, transcripts, by='target_id', all.x=TRUE)
tximp_tx_results <- mutate(tximp_tx_results, genes = gene_id)
tximp_tx_results <- select(tximp_tx_results, genes, pvalue, baseMean)
tximp_tx_results <- tximp_tx_results %>% group_by(genes) %>%
	summarise(tx_deseq = lancaster(pvalue, baseMean),
		tx_deseq_f = fisher(pvalue),
		tx_deseq_min = MinMethod(pvalue))
summary_table <- merge(summary_table, tximp_tx_results, by='genes', all.x=TRUE)

#limma on transcripts
print('limma on tx')
limma_results <- limma(tx_counts, design)
limma_summary <- data.frame(pvalue = limma_results$p.value[,2], mean = limma_results$Amean, target_id = rownames(tx_counts))
limma_summary <- merge(limma_summary, transcripts, by='target_id', all.x=TRUE)
limma_summary <- mutate(limma_summary, genes = gene_id)
limma_summary <- limma_summary %>% group_by(genes) %>%
	summarise(tx_limma = lancaster(pvalue, exp(mean)),
		tx_limma_f = fisher(pvalue),
		tx_limma_min = MinMethod(pvalue))
summary_table <- merge(summary_table, limma_summary, by='genes', all.x=TRUE)

#sleuth pipeline on transcript counts
s2c$path <- files
so <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE) 
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_wt(so, 'condition1', 'full')
sleuth_table <- sleuth_results(so, 'condition1', which_model = 'full')
if(nrow(sleuth_table)!= nrow(transcripts))
{
	print('error, number of rows in sleuth table and index does not match')
}
sleuth_table <- merge(sleuth_table, transcripts, by='target_id', all.x=TRUE)
sleuth_summary <- sleuth_table %>% group_by(gene_id) %>%
	summarise(tx_sleuth = lancaster(pval, exp(mean_obs)),
		tx_sleuth_f = fisher(pval),
		tx_sleuth_min = MinMethod(pval))
sleuth_summary <- mutate(sleuth_summary, genes = gene_id)
sleuth_summary <- select(sleuth_summary, genes, tx_sleuth, tx_sleuth_f, tx_sleuth_min)
if(nrow(summary_table) != nrow(sleuth_summary))
{
	print('error, number of rows should be equal in summary table and sleuth summary')
}
summary_table <- merge(summary_table, sleuth_summary, by='genes', all.x=TRUE)

print('finished transcript analysis')
print(dim(summary_table))
saveRDS(summary_table, file.path(experiment, exp_string, 'summary_table2.rds'))

sims_info <- readRDS(file.path(experiment2, 'sims.rds'))
sims_info <- sims_info[[exp_num]]$info
sims_info <- merge(sims_info, transcripts, by='target_id', left.all = TRUE)
sims_info  <- mutate(sims_info, de = is_de, genes = gene_id)
sims_info <- select(sims_info, genes, de)
sims_info <- sims_info %>% group_by(genes) %>% summarise(de = any(de))
summary_table <- merge(summary_table, sims_info, by='genes', all.x = TRUE)

rds <- file.path(experiment, exp_string, 'summary_table.rds')
saveRDS(summary_table, rds)

