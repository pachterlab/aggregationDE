
source('roc_curve.R')
args <- commandArgs(TRUE)
#base_path <- '~/geuvadis/gfr_3_3_20_42_2/'
#base_path <- '~/geuvadis/isoform_3_3_20_1_1/'
#base_path <- '~/geuvadis/gcd_3_3_20_1_2/'
base_path <- args[[1]]

paths <- sapply(1:20, function(i) file.path(base_path, 'lrt', paste0('exp_',i), 'summary_table.rds'))

col_names <- c('gene_sleuth', 'tx_sleuth', 'tx_sleuth_min', 'tcc_sleuth')
titles <- c('sleuth - Gene', 'sleuth - Lancaster Tx', 'sleuth - Sidak Tx', 'sleuth - TCC')

#select and adjust columns
summary_table_list <- lapply(paths, readRDS)
selected <- lapply(summary_table_list, function(x) dplyr::select(x, col_names))
adjusted <- lapply(selected, function(x)
				 {
				 		 apply(x, 2, function(y) p.adjust(y, 'BH'))
				 		 x$genes <- summary_table_list[[1]]$genes
				 })
#calculate fdr/sensitivity per sim and average
average_fdrs <- lapply(seq_along(col_names), function(i)
	averaging(adjusted, col_names[i], titles[i]))
if(length(col_names) != length(average_fdrs))
{
	print('error, lengths are not equal')
}
saveRDS(average_fdrs, file.path(base_path, 'lrt', 'average_fdrs.rds'))

