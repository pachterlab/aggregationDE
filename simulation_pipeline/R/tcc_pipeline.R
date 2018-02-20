library(dplyr)
library(DESeq2)
library(sleuth)
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

# need ncol counts == nrow design == nsamples
#print('starting DESeq on TCC analysis')
DESeqData <- DESeq_lrt(t(df), design)
deseq_table <- as.data.frame(results(DESeqData))
saveRDS(deseq_table, file.path(experiment, 'lrt', exp_string, 'deseq_tcc_table.rds'))
deseq_pvalues <- deseq_table$pvalue

#sleuth_tcc
print('starting sleuth on TCC analysis')
tcc_files <- sapply(1:6, function(i) file.path(experiment, exp_string, 'sleuth_tcc', i, 'abundance.h5'))
print(tcc_files)
sample <- c('s1','s2','s3','s4','s5','s6')
s2c <- data.frame(sample, condition = design, path = tcc_files, stringsAsFactors = FALSE)
print(s2c)
so <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE, filter_fun = filter_func)
sleuth_table <- sleuth_lrt_pipeline(so) 
sleuth_table$target_id <- as.numeric(sleuth_table$target_id)
sleuth_table <- sleuth_table[order(sleuth_table$target_id),] 
sleuth_pvalues <- sleuth_table$pval
saveRDS(sleuth_table, file.path(experiment, 'lrt', exp_string, 'sleuth_tcc_table.rds'))

#map tccs to genes
tcc2genemap <- make_tcc_to_gene_map2(transcripts, ec_map)
tcc_table <- data.frame(genes = tcc2genemap, dpvalue = deseq_pvalues, spvalue = sleuth_pvalues, meancounts = colMeans(df))
saveRDS(tcc_table, file.path(experiment, 'lrt', exp_string,'tcc_table.rds'))

summary_table <- tcc_table %>% group_by(genes) %>%
	summarise(n = n(),
		tcc_deseq = lancaster(dpvalue, normalize_tccs(meancounts)),
		tcc_sleuth = lancaster(spvalue, normalize_tccs(meancounts)))
summary_table <- summary_table[-1,]

print('finished TCC analyses')
print(head(summary_table))
saveRDS(summary_table, file.path(experiment, 'lrt', exp_string, 'tcc_summary_table.rds'))
