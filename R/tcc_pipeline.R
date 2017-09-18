
library(sleuth)
source('aggregation.R')
#change path to where your reads are stored
base_path <- '/home/lynnyi/SRP091911'

print('reading ec map')
#matrix.ec is product of pseudoquant
ec_map_path <- file.path(base_path, 'pseudoquant/matrix.ec')
print(ec_map_path)
#ec_map is product of making kallisto index that maps EC number to transcript number
ec_map <- read.table(ec_map_path, sep='\t', stringsAsFactors=FALSE)
colnames(ec_map) <- c('tcc_num', 'transcript_num')

print('reading matrix')
#reading matrix of pseudoquant, an output of pseudoquant
in_matrix <- file.path(base_path, 'pseudoquant/matrix.tsv')
print(in_matrix)
df <- scan(in_matrix, skip = 1, what=numeric(), sep='\t')
#nrow = number of your samples in matrix 
df <- matrix(df, nrow=21)
df <- df[-1,]
print(dim(df))

#change your gene to transcript mapping
print('reading transcript_gene map')
transcripts <- read.table('/home/lynnyi/transcriptomes/Rattus_norvegcius.Rnor_6.0.transcripts', stringsAsFactors=FALSE)
colnames(transcripts) <- c('target_id', 'genes')

#removing everything after .
transcripts$target_id <- gsub('\\..*', '', transcripts$target_id)
transcripts$genes <- gsub('\\..*', '', transcripts$genes)

#read sample table
s2c <- read.table(file.path(base_path,'/simple_sample_table.txt'), sep = '\t', header=TRUE)
names <- s2c$Run_s
paths <- file.path(base_path, '/sleuth_tcc/', names, 'abundance.h5')
s2c$sample <- names
s2c$path <- paths

filter_func <- function(row)
{
		    mean(row) > 5
}

#run sleuth pipeline and change variable you're testing
so <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE, filter_fun = filter_func)
so <- sleuth_fit(so, ~stretch, 'full')
so <- sleuth_wt(so, 'stretchstretch', 'full')

sleuth_table <- sleuth_results(so, 'stretchstretch', which_model = 'full')
sleuth_table$target_id <- as.numeric(sleuth_table$target_id)
sleuth_table <- sleuth_table[order(sleuth_table$target_id),]

print('finished sleuth pipeline')
head(sleuth_table)
saveRDS(sleuth_table, file.path(base_path, 'sleuth_tcc_table.rds'))

#tcc2genemap takes ec_map and provides list of unambiguous gene names
print('starting maps')
tcc2genemap <- make_tcc_to_gene_map(transcripts, ec_map)
tcc_table <- data.frame(genes = tcc2genemap, sleuth_pval = sleuth_table$pval, meancounts = colMeans(df))
saveRDS(tcc_table, file.path(base_path, 'tcc_table.rds'))

#doing gene level DE by aggregating pvalues
gene_de_table <- tcc_table %>% group_by(genes) %>%
		summarise(pval = lancaster(sleuth_pval, meancounts))
gene_de_table <- gene_de_table[-1,]
saveRDS(gene_de_table,file.path(base_path,'tcc_pipeline_results.rds'))

