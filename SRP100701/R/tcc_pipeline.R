library(sleuth)
library(aggregation)
source('~/TCC/misc.R')

base_path <- '/home/lynnyi/SRP100701'
print('reading ec map')
ec_map_path <- file.path(base_path, 'pseudoquant/matrix.ec')
ec_map <- read.table(ec_map_path, sep='\t', stringsAsFactors=FALSE)
colnames(ec_map) <- c('tcc_num', 'transcript_num')

print('reading matrix')
in_matrix <- file.path(base_path, 'pseudoquant/matrix.tsv')
print(in_matrix)
df <- scan(in_matrix, skip = 1, what=numeric(), sep='\t')
df <- matrix(df, nrow=25)
df <- df[-1,]
print(dim(df))

print('reading transcript_gene map')
transcripts <- read.table('/home/lynnyi/transcriptomes/Mus_musculus.GRCm38.cdna.rel.88.transcripts', stringsAsFactors=FALSE)
colnames(transcripts) <- c('target_id', 'genes')
transcripts$target_id <- gsub('\\..*', '', transcripts$target_id)
transcripts$genes <- gsub('\\..*', '', transcripts$genes)

s2c <- read.table(file.path(base_path,'/sample_table.txt'), sep = '\t', header=TRUE)
names <- s2c$Run_s
paths <- file.path(base_path, '/sleuth_tcc/', names, 'abundance.h5')
s2c$sample <- names
s2c$path <- paths
s2c$gender <- s2c$gender_s
s2c$treatment <- s2c$treatment_s
s2c$region <- s2c$tissue_region_s


so <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE, filter_fun = filter_func)
so <- sleuth_fit(so, ~gender+treatment+region, 'full')
so <- sleuth_fit(so, ~gender+region, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table <- sleuth_results(so, 'reduced:full', test_type='lrt')

sleuth_table$target_id <- as.numeric(sleuth_table$target_id)
sleuth_table <- sleuth_table[order(sleuth_table$target_id),]

print('finished sleuth pipeline')
head(sleuth_table)
saveRDS(sleuth_table, file.path(base_path, 'wt', 'sleuth_table_tcc.rds'))

print('starting maps')
ptm <- proc.time()
tcc2genemap <- make_tcc_to_gene_map2(transcripts, ec_map)
proc.time() - ptm

tcc_table <- data.frame(genes = tcc2genemap, sleuth_pval = sleuth_table$pval, meancounts = colMeans(df))
saveRDS(tcc_table, file.path(base_path, 'wt','tcc_table.rds'))

gene_de_table <- tcc_table %>% group_by(genes) %>%
		summarise(pval = lancaster(sleuth_pval, normalize_tccs(meancounts)))
gene_de_table <- gene_de_table[-1,]
gene_de_table$qval <- p.adjust(gene_de_table$pval, 'BH')
saveRDS(gene_de_table,file.path(base_path,'wt','tcc_results.rds'))

