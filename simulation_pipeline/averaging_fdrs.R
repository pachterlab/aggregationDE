
source('roc_curve.R')


args <- commandArgs(TRUE)

#base_path <- '~/geuvadis/gfr_3_3_20_42_2/'
#base_path <- '~/geuvadis/isoform_3_3_20_1_1/'
#base_path <- '~/geuvadis/gcd_3_3_20_1_2/'
base_path <- args[[1]]

paths <- sapply(1:20, function(i) file.path(base_path, paste0('exp_',i), 'summary_table.rds'))


col_names <- c('gene_sleuth', 'tx_sleuth', 'tx_sleuth_min', 'tcc_sleuth')
titles <- c('Sleuth - Gene', 'Sleuth - Lancaster Tx', 'Sleuth - Min Tx', 'Sleuth - TCC')

summary_table_list <- lapply(paths, readRDS)

adjust <- function(summary_table_list)
{
	for(i in 1:length(summary_table_list))
	{
		for(j in 3:17)
		{
			print(i)
			print(j)
			summary_table_list[[i]][,j] <-
				p.adjust(summary_table_list[[i]][,j], method='BH')
		}
	}
	summary_table_list
}

summary_table_list <- adjust(summary_table_list)

saveRDS(summary_table_list, file.path(base_path, 'adjusted_summary_table.rds'))

average_fdrs <- lapply(seq_along(col_names), function(i)
	averaging(summary_table_list, col_names[i], titles[i]))

if(length(col_names) != length(average_fdrs))
{
	print('error, lengths are not equal')
}

saveRDS(average_fdrs, file.path(base_path, 'average_fdrs.rds'))


plot_fdrs(average_fdrs[[1]], average_fdrs[[2]], average_fdrs[[3]], average_fdrs[[4]], file = file.path(base_path, 'average_fdrs.png'), title =  base_path)


