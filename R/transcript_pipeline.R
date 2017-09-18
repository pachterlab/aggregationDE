
library(sleuth)
source('../R/aggregation.R')

#change path and transcript to gene mapping table
base_path <- '/home/lynnyi/SRP091911'

print('reading transcript_gene map')
transcripts <- read.table('/home/lynnyi/transcriptomes/Rattus_norvegcius.Rnor_6.0.transcripts', stringsAsFactors=FALSE)
colnames(transcripts) <- c('target_id', 'genes')
#transcripts$target_id <- gsub('\\..*', '', transcripts$target_id)
transcripts$genes <- gsub('\\..*', '', transcripts$genes)

s2c <- read.table(file.path(base_path,'/simple_sample_table.txt'), sep = '\t', header=TRUE, stringsAsFactors=FALSE)
names <- s2c$Run_s
paths <- file.path(base_path, '/kallisto/', names, 'abundance.h5')
s2c$sample <- names
s2c$path <- paths

so <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE)
so <- sleuth_fit(so, ~stretch, 'full')
so <- sleuth_wt(so, 'stretchstretch', 'full')
sleuth_table <- sleuth_results(so, 'stretchstretch', which_model = 'full')
print('finished sleuth pipeline')
head(sleuth_table)

saveRDS(sleuth_table, file.path(base_path, 'sleuth_transcript_table.rds'))

merged_results <- merge(sleuth_table, transcripts, by='target_id', all.x=TRUE)
results <- results %>% group_by(genes) %>% summarise(lan = lancaster(pval, exp(mean_obs)), weight = sum(exp(mean_obs), na.rm=TRUE), min = MinMethod(pval))

results$lan.adjust <- p.adjust(results$lan, method = 'BH')
results$min.adjust <- p.adjust(results$min, method = 'BH')

saveRDS(results, file.path(base_path, 'transcript_pipeline_results.rds'))

