
library(sleuth)
source('../R/aggregation.R')

#change path to where your reads are stored
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

print('first 10 rows')
df <- df[c(1:10),]


print('reading transcript_gene map')
#change your gene to transcript mapping
transcripts <- read.table('/home/lynnyi/transcriptomes/Rattus_norvegcius.Rnor_6.0.transcripts', stringsAsFactors=FALSE)
colnames(transcripts) <- c('target_id', 'genes')
transcripts$target_id <- gsub('\\..*', '', transcripts$target_id)
transcripts$genes <- gsub('\\..*', '', transcripts$genes)

s2c <- read.table(file.path(base_path,'/simple_sample_table.txt'), sep = '\t', header=TRUE)
names <- s2c$Run_s
paths <- file.path(base_path, '/sleuth_tcc/', names, 'abundance.h5')
s2c$sample <- names
s2c$path <- paths

filter_func <- function(row)
{
		    mean(row) > 5
}

so <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE, filter_fun = filter_func)
so <- sleuth_fit(so, ~stretch, 'full')
so <- sleuth_wt(so, 'stretchstretch', 'full')

sleuth_table <- sleuth_results(so, 'stretchstretch', which_model = 'full')
sleuth_table$target_id <- as.numeric(sleuth_table$target_id)
sleuth_table <- sleuth_table[order(sleuth_table$target_id),]

print('finished sleuth pipeline')
head(sleuth_table)
saveRDS(sleuth_table, file.path(base_path, 'sleuth_tcc_table.rds'))

print('starting maps')
tcc2genemap <- make_tcc_to_gene_map(transcripts, ec_map)

tcc_table <- data.frame(genes = tcc2genemap, sleuth_pval = sleuth_table$pval, meancounts = colMeans(df))
saveRDS(tcc_table, file.path(base_path, 'tcc_table.rds'))

#doing gene level DE
gene_de_table <- tcc_table %>% group_by(genes) %>%
		summarise(pval = lancaster(sleuth_pval, meancounts))
gene_de_table <- gene_de_table[-1,]
saveRDS(gene_de_table,file.path(base_path,'tcc_pipeline_results.rds'))

