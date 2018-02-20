
#plot p-values to generate Sup Fig 4

base_path = '~/SRP100701/'
sleuth_table <- readRDS(file.path(base_path, '/lrt/', 'sleuth_table_tx.rds'))

png(file.path(base_path, 'lrt', 'plots', 'uniformity.png'), height = 2500, width = 2500, res = 1000)

hist(sleuth_table$pval, breaks=20, xlim=c(0,1))

dev.off()
