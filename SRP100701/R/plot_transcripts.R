library(dplyr)
library(reshape)
library(sleuth)
base_path = '~/SRP100701/'
gene_so <- readRDS(file.path(base_path, 'wt', 'gene_so.rds'))
gene_counts <- gene_so$obs_norm
tx_so <- readRDS(file.path(base_path, 'wt', 'tx_so.rds'))
tx_counts <- tx_so$obs_norm
sleuth_table <- readRDS(file.path(base_path, 'wt', 'tx_sleuth_table.rds'))

s2c <- read.table(file.path(base_path,'/sample_table.txt'), sep = '\t', header=TRUE, stringsAsFactors=FALSE)
s2c$sample <- s2c$Run_s
s2c$treatment <- s2c$treatment_s
condition1 <- filter(s2c, treatment == 'Vehicle')
condition1 <- condition1$sample
condition2 <- filter(s2c, treatment == 'Dexamethasone')
condition2 <- condition2$sample

t2g <- read.table('/home/lynnyi/transcriptomes/Mus_musculus.GRCm38.cdna.rel.88.transcripts', sep='\t', stringsAsFactors=FALSE)
colnames(t2g) <- c('target_id', 'genes')

gene_pipeline_results <- readRDS(file.path(base_path, 'wt', 'gene_results.rds'))
tx_pipeline_results <- readRDS(file.path(base_path, 'wt', 'tx_results.rds'))

#genes_to_plot
gene_pipeline_results <- select(gene_pipeline_results, target_id, pval, qval)
results <- merge(gene_pipeline_results, tx_pipeline_results, by.x='target_id', by.y ='genes', all=TRUE)
#interesting_results <- filter(results, lan < .05, pval > .05, min  < .05)
#interesting_results <- interesting_results[order(interesting_results$lan),]
#interesting_genes <- interesting_results$target_id

get_tpm_obs <- function(obs_counts, gene_counts, sleuth_table, transcripts, condition1, condition2, gene, gene_name, gene_pval, lancaster)
{

	x <- filter(gene_counts, target_id == gene, sample %in% condition1)
	y <- filter(gene_counts, target_id == gene, sample %in% condition2)
	m <- c(mean(x$tpm), mean(y$tpm))
	se <- c(sd(x$tpm)/sqrt(length(x$tpm)), sd(y$tpm)/sqrt(length(y$tpm)))
	m <- melt(m)
	print(m)
	print(m)
	se <- melt(se)
	print(se)
	upper <- m + se
	lower <- m - se
	treatment <- c('Vehicle', 'Dexamethasone')
	transcript <- c(gene)
	m  <- cbind(treatment, transcript, m, se, upper, lower)
	colnames(m) <- c('Treatment', 'Transcript', 'mean', 'se', 'upper', 'lower')
	m$pval <- gene_pval
	print(m)
	
	table <- sapply(transcripts, function(transcript)
		   { x <- filter(obs_counts, (target_id) ==  transcript, (sample) %in% condition1)
			 y <- filter(obs_counts, (target_id) == transcript, sample %in% condition2)
		     counts_mean <- c( mean(x$tpm), mean(y$tpm))
		     counts_se <- c(sd(x$tpm)/sqrt(length(x$tpm)), sd(y$tpm)/sqrt(length(x$tpm)))
			 names(counts_mean) <- c('Vehicle', 'Dexamethasone')
			 names(counts_se) <- c('Vehicle', 'Dexamethasone')
			 c(counts_mean, counts_se)
		   })
	
	counts <- table[1:2,]
	se <- table[3:4,] 
	counts <- melt(counts)
	se <- melt(se)
	colnames(counts) <- c('Treatment', 'Transcript', 'mean')
	colnames(se) <- c('Treatment', 'Transcript', 'se')
	counts <- merge(counts, se, by=c('Treatment','Transcript'))
	
	counts$upper <- counts$mean + counts$se
	counts$lower <- counts$mean - counts$se

	sleuth_table <- select(sleuth_table, target_id, pval)
	counts <- merge(counts, sleuth_table, by.x = 'Transcript', by.y='target_id')
	counts$pval <- signif(counts$pval, 2)

	##change this
	#counts <-rbind(counts, m)
	counts <- m	
	counts$Transcript <- factor(counts$Transcript, levels(counts$Transcript))
	counts$Treatment <- factor(counts$Treatment, levels(counts$Treatment)[c(2,1)])
	
	paste_new <- function(x, y) paste(x, y, sep=', ')
	new_labels <- mapply(paste_new, counts$Transcript, counts$pval)
	new_labels <- unique(new_labels)
	print(new_labels)	

	#default coloring/labeling
	p <- ggplot(counts, aes(x=Transcript , y=mean, fill=Treatment)) + geom_bar(stat = "identity", position = 'dodge') + ylab('Mean TPM') + labs(fill = gene) + xlab(lancaster)
	p <- p + scale_x_discrete(labels = unique(new_labels))
	p <- p + geom_errorbar(aes(ymax = counts$upper, ymin=counts$lower), position = 'dodge')
	p <- p + theme(text = element_text(size=11)) 
	png(paste0(base_path, '/wt/plots/', gene, '.png'), height = 5000, width = 5000, res = 1000)
	print('executing')
	print(p)
	dev.off()
}


get_obs <- function(obs_counts, gene_counts, sleuth_table, transcripts, condition1, condition2, gene, gene_name, gene_pval, lancaster)
{
	#get gene count values
	x <- filter(gene_counts, target_id == gene, sample %in% condition1)
	y <- filter(gene_counts, target_id == gene, sample %in% condition2)
	m <- c(mean(x$scaled_reads_per_base), mean(y$scaled_reads_per_base))
	print(x$scaled_reads_per_base)
	print(x)
	print(var(x$scaled_reads_per_base))
	print(var(y$scaled))
	se <- c(sd(x$scaled)/sqrt(length(x$scaled)), sd(y$scaled)/sqrt(length(y$scaled)))
	m <- melt(m)
	se <- melt(se)
	print(se)
	upper <- m + se
	lower <- m - se
	treatment <- c('Vehicle', 'Dexamethasone')
	transcript <- c(gene)
	m  <- cbind(treatment, transcript, m, se, upper, lower)
	colnames(m) <- c('Treatment', 'Transcript', 'mean', 'se', 'upper', 'lower')
	m$pval <- gene_pval
	print(m)

	#get transcript count values
	table <- sapply(transcripts, function(transcript)
		   { x <- filter(obs_counts, (target_id) ==  transcript, (sample) %in% condition1)
			 y <- filter(obs_counts, (target_id) == transcript, sample %in% condition2)
		     counts_mean <- c( mean(x$est_counts), mean(y$est_counts))
			 print(var(x$est_counts))
			 print(var(y$est_counts))
			 counts_se <- c(sd(x$est_counts)/sqrt(length(x$est_counts)), sd(y$est_counts)/sqrt(length(x$est_counts)))
			 names(counts_mean) <- c('Vehicle', 'Dexamethasone')
			 names(counts_se) <- c('Vehicle', 'Dexamethasone')
			 c(counts_mean, counts_se)
		   })
	
	counts <- table[1:2,]
	se <- table[3:4,] 
	counts <- melt(counts)
	se <- melt(se)
	colnames(counts) <- c('Treatment', 'Transcript', 'mean')
	colnames(se) <- c('Treatment', 'Transcript', 'se')
	counts <- merge(counts, se, by=c('Treatment','Transcript'))
	
	counts$upper <- counts$mean + counts$se
	counts$lower <- counts$mean - counts$se

	sleuth_table <- select(sleuth_table, target_id, pval)
	counts <- merge(counts, sleuth_table, by.x = 'Transcript', by.y='target_id')
	counts$pval <- signif(counts$pval, 2)

	#bind gene and transcript count values
	counts <-rbind(counts, m)
	#counts <- m
	counts$Transcript <- factor(counts$Transcript, levels(counts$Transcript))
	counts$Treatment <- factor(counts$Treatment, levels(counts$Treatment)[c(2,1)])
	#sleuth table
	paste_new <- function(x, y) paste(x, y, sep=', ')
	new_labels <- mapply(paste_new, counts$Transcript, counts$pval)
	new_labels <- unique(new_labels)
	print(new_labels)	

	#default coloring/labeling
	p <- ggplot(counts, aes(x=Transcript , y=mean, fill=Treatment)) + geom_bar(stat = "identity", position = 'dodge') + ylab('Mean Counts') + labs(fill = gene) + xlab(lancaster)
	p <- p + scale_x_discrete(labels = unique(new_labels))
	p <- p + geom_errorbar(aes(ymax = counts$upper, ymin=counts$lower), position = 'dodge')
	p <- p + theme(text = element_text(size=11))
	png(paste0(base_path, '/wt/plots/', gene, '.png'), height = 5000, width = 5000, res = 1000)
	print(p)
	dev.off()
}

#sleuth_table is transcript table
#condition1 and condition2 refers to samples
plot_me <- function(interesting_genes) {
  sapply(interesting_genes, function(gene)
	   {
	   		   transcripts <- filter(t2g, genes == gene)$target_id
	   		   st <- filter(sleuth_table, target_id %in% transcripts)
	   		   gene_pval <- filter(gene_pipeline_results, target_id == gene)$pval
	   		   gene_pval <- signif(gene_pval, 2)
	   		   lan_pval <- filter(tx_pipeline_results, genes == gene)$lan
	   		   lan_pval <- signif(lan_pval, 2)
	   		   label <- paste0('Lancaster p-val: ', lan_pval)
	   		   gene_name <- paste(gene, gene_pval, label, sep=', ')
	   		   get_obs(tx_counts, gene_counts, st, transcripts, condition1, condition2, gene, gene_name, gene_pval, label)

	   })
}

interesting_genes <- c('ENSMUSG00000001305.5', 'ENSMUSG00000016409.4', 'ENSMUSG00000022456.15', 'ENSMUSG00000062190.12')
plot_me(interesting_genes)
