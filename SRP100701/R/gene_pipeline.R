
#gene level sleuth analysis
library(sleuth)
base_path <- '/home/lynnyi/SRP100701'

t2g <- read.table('/home/lynnyi/transcriptomes/Mus_musculus.GRCm38.cdna.rel.88.transcripts', stringsAsFactors=FALSE)
colnames(t2g) <- c('target_id', 'genes')

s2c <- read.table(file.path(base_path,'/sample_table.txt'), sep = '\t', header=TRUE)
names <- s2c$Run_s
paths <- file.path(base_path, '/kallisto/', names, 'abundance.h5')
s2c$sample <- names
s2c$path <- paths
s2c$gender <- s2c$gender_s
s2c$treatment <- s2c$treatment_s
s2c$region <- s2c$tissue_region_s

so <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE, target_mapping = t2g, aggregation_column = 'genes')
print('finished importing data')
so <- sleuth_fit(so, ~gender+treatment+region, 'full')
#so <- sleuth_fit(so, ~gender+region, 'reduced')
#so <- sleuth_wt(so, 'reduced', 'full')
#sleuth_table <- sleuth_results(so, 'reduced:full', 'wt')
so <- sleuth_wt(so, 'treatmentVehicle', 'full')
sleuth_table <- sleuth_results(so, 'treatmentVehicle', which_model='full')
print('finished sleuth pipeline')


saveRDS(sleuth_table, file.path(base_path, 'wt', 'gene_results.rds'))
saveRDS(so, file.path(base_path, 'wt', 'gene_so.rds'))
