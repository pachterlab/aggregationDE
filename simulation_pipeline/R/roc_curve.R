
library(dplyr)

calculate_fdr <- function(true_DE, qvalues, title)
{
	df <- data.frame(true_DE, qvalues)
	df <- df[order(df$qvalues), ]
	total_positive <- sum(df$true_DE)
	total_negative <- nrow(df) - total_positive
	df <- filter(df, !is.na(df$qvalues))
	print(head(df))
	fdr <- sapply(seq(df$qvalues), function(i) sum(!df$true_DE[1:i]) / i)
	sensitivity <- sapply(seq(df$qvalues), function(i) sum(df$true_DE[1:i])/total_positive)
	n_found <- sapply(seq(df$qvalues), function(i) sum(df$true_DE[1:i]))
	n_de <- (1:nrow(df))
	five <- min(which(df$qvalues > .05)) -1 
	ten <- min(which(df$qvalues > .1)) - 1
	twenty <- min(which(df$qvalues > .2)) -1
	
	list(fdr = fdr, sensitivity = sensitivity,
		n_found = n_found,
		n_de = n_de,
		total_positive = total_positive,
		total_negative = total_negative, 
		qvalues = df$qvalues,
		five=five, ten=ten, twenty=twenty, title=title)
}


average_fdr <- function(fdr_list)
{
	#convert list to data frame
	fdrs <- lapply(fdr_list, function(x)
		data.frame(fdr = x$fdr,
		sensitivity = x$sensitivity,
		n_de = x$n_de,
		n_found = x$n_found,
		qvalues = x$qvalues,
		total_positive = x$total_positive,
		total_negative = x$total_negative))
	fdrs <- do.call(rbind, fdrs)
	print(head(fdrs))
	print(dim(fdrs))
	#perform averaging
	average_fdr <- fdrs %>% group_by(n_de) %>%
		summarise(fdr = mean(fdr),
			qvalues = mean(qvalues),
			sd_sensitivity= sd(sensitivity),
			sensitivity = mean(sensitivity),
			n_found = mean(n_found),
			total_positive = mean(total_positive),
			total_negative = mean(total_negative))
	print(dim(average_fdr))
	print(head(average_fdr))	
	five <- min(which(average_fdr$qvalues > .05)) -1
    ten <- min(which(average_fdr$qvalues > .1)) - 1
    twenty <- min(which(average_fdr$qvalues > .2)) -1
	
	list(fdr = average_fdr$fdr,
		sensitivity = average_fdr$sensitivity,
		sd_sensitivity = average_fdr$sd_sensitivity,
		n_found = average_fdr$n_found,
		n_de = average_fdr$n_de,
		qvalues = average_fdr$qvalues,
		five = five, ten=ten, twenty=twenty,
		title = fdr_list[[1]]$title,
		total_positive = average_fdr$total_positive,
		total_negative = average_fdr$total_negative)
}

plot_fdr <- function(fdr, color)
{
	lines(fdr$fdr, fdr$sensitivity, 'l', col=color)
	idx <- c(fdr$five, fdr$ten, fdr$twenty)
	points(fdr$fdr[idx], fdr$sensitivity[idx], pch = c(1,2,0), cex=1)
}


plot_fdrs <- function(fdr1, fdr2 = NULL, fdr3=NULL, fdr4=NULL, fdr5 = NULL,
	fdr6 =NULL, fdr7 = NULL, fdr8 = NULL, file, title)
{
	png(file)
	plot(fdr1$fdr, fdr1$sensitivity, col = 'dodgerblue3', type='l')
	idx <- c(fdr1$five, fdr1$ten, fdr1$twenty)
        points(fdr1$fdr[idx], fdr1$sensitivity[idx], pch = c(1,2,0), cex=1)
	if(!is.null(fdr2))
	{
		plot_fdr(fdr2, 'indianred3')
	}
	if(!is.null(fdr3))
	{
		plot_fdr(fdr3, 'orchid4')
	}
	if(!is.null(fdr4))
	{
		plot_fdr(fdr4, 'mediumseagreen')
	}
	if(!is.null(fdr5))
	{
		plot_fdr(fdr5, 'orange1')
	}
	if(!is.null(fdr6))
	{
		plot_fdr(fdr6, 'plum')
	}
	if(!is.null(fdr7))
	{
		plot_fdr(fdr7, 'steelblue2')
	}
	if(!is.null(fdr8))
	{
		plot_fdr(fdr8, 'purple')
	}
	title(title)
	
	legend('bottomright', c(fdr1$title, fdr2$title, fdr3$title, fdr4$title, fdr5$title, fdr6$title, fdr7$title, fdr8$title), lty=rep(1,8), lwd=rep(2.5, 8), col=c('dodgerblue4', 'indianred3', 'orchid4', 'mediumseagreen', 'orange1', 'plum', 'steelblue2', 'purple'))
	dev.off()
}

averaging <- function(summary_table_list, colname, title)
{
	i <- which(colnames(summary_table_list[[1]]) == colname)
	print('column')
	print(i)
	print('calculating fdrs for each summary_table')
	fdrs_list <- lapply(summary_table_list, function(x)
		calculate_fdr(true_DE = x$de,
				qvalues = x[,i],
				title = title))
	print('running average_fdr on fdr_list')
	average <- average_fdr(fdrs_list)
	average
}

