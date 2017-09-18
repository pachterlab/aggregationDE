
source('../mamabear/R/de_benchmark.R')
library(ggplot2)

#base_path <- '~/geuvadis/gfr_3_3_20_42_2'
#base_path <- '~/geuvadis/isoform_3_3_20_1_1/'
#base_path <- '~/geuvadis/gcd_3_3_20_1_2/'

args <- commandArgs(TRUE)
base_path <- args[[1]]

convert_to_bear <- function(fdr)
{
	bear_fdr <- data.frame(
		sensitivity = fdr$sensitivity,
		true_fdr = fdr$fdr,
		#pvalues already adjusted
		qval = fdr$pvalues,
		method = fdr$title,
		nde = fdr$n_de,
		sd_sensitivity = fdr$sd_sensitivity,
		p = fdr$total_positive, 
		n = fdr$total_negative
		)
}

fdrs <- readRDS(file.path(base_path, 'average_fdrs.rds'))
bear_fdrs <- lapply(fdrs, convert_to_bear)
bear_fdrs <- do.call(rbind, bear_fdrs)

p <- fdr_efdr_power_plot(bear_fdrs, isolines = FALSE)

png(file.path('~/mama_plot_fa.png'), height = 8000, width = 8000, res=1200)
p <- p + coord_cartesian(xlim = c(-0.05, 1), ylim = c(-0.05, 1), expand = FALSE)
p <- p + theme_bw()

print(p)
dev.off()

png(file.path('~/mama_plot_fb.png'), height = 8000, width=8000, res=1200)
p <- p + coord_cartesian(xlim=c(-.01, .15), ylim=c(-.01,.4))
p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(p)
dev.off()


