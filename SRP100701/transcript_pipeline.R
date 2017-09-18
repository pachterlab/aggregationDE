
library(sleuth)
source('../R/aggregation.R')

#change path to where data is, change transcript to gene mapping location
base_path <- '/home/lynnyi/SRP100701'
s2c <- read.table(file.path(base_path,'/sample_table.txt'), sep = '\t', header=TRUE, stringsAsFactors=FALSE)
names <- s2c$Run_s
paths <- file.path(base_path, '/kallisto/', names, 'abundance.h5')
s2c$sample <- names
s2c$path <- paths
s2c$gender <- s2c$gender_s
s2c$treatment <- s2c$treatment_s
s2c$region <- s2c$tissue_region_s

so <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE)
so <- sleuth_fit(so, ~gender+treatment+region, 'full')
so <- sleuth_wt(so, 'treatmentVehicle', 'full')
sleuth_table <- sleuth_results(so, 'treatmentVehicle', which_model = 'full')
print('finished sleuth pipeline')
head(sleuth_table)

saveRDS(sleuth_table, file.path(base_path, 'sleuth_transcript_table.rds'))

transcripts <- read.table('/home/lynnyi/transcriptomes/Mus_musculus.GRCm38.cdna.rel.88.transcripts', sep='\t', stringsAsFactors=FALSE)
colnames(transcripts) <- c('target_id', 'genes')

results <- merge(sleuth_table, transcripts, by='target_id', all.x=TRUE)
results <- results %>% group_by(genes) %>% summarise(lan = lancaster(pval, exp(mean_obs)), weight = sum(exp(mean_obs), na.rm=TRUE), min = MinMethod(pval))

results$lan.adjust <- p.adjust(results$lan, method = 'BH')
results$min.adjust <- p.adjust(results$min, method = 'BH')

saveRDS(results, file.path(base_path, 'transcript_pipeline_results.rds'))

