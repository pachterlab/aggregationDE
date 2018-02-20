library(sleuth)
library(aggregation)
source('~/TCC/misc.R')

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
#so <- sleuth_fit(so, ~gender+region, 'reduced')
#so <- sleuth_wt(so, 'reduced', 'full')
#sleuth_table <- sleuth_results(so, 'reduced:full', test_type='wt')
so <- sleuth_wt(so, 'treatmentVehicle', 'full')
sleuth_table <- sleuth_results(so, 'treatmentVehicle', which_model = 'full')

print('finished sleuth pipeline')
head(sleuth_table)

saveRDS(sleuth_table, file.path(base_path, 'wt', 'tx_sleuth_table.rds'))
transcripts <- read.table('/home/lynnyi/transcriptomes/Mus_musculus.GRCm38.cdna.rel.88.transcripts', sep='\t', stringsAsFactors=FALSE)
colnames(transcripts) <- c('target_id', 'genes')
results <- merge(sleuth_table, transcripts, by='target_id', all.x=TRUE)
results <- results %>% group_by(genes) %>% summarise(lan = lancaster(pval, mean_obs), weight = sum(mean_obs, na.rm=TRUE), min = sidak(pval))

results$lan.adjust <- p.adjust(results$lan, method = 'BH')
results$min.adjust <- p.adjust(results$min, method = 'BH')
saveRDS(so, file.path(base_path, 'wt', 'tx_so.rds'))
saveRDS(results, file.path(base_path, 'wt', 'tx_results.rds'))
