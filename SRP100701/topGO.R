
library(topGO)
library(stats)
library(dplyr)
library(ggplot2)
source('../R/aggregation.R')
base_dir <- '/home/lynnyi/SRP100701'

print('reading files')
results <- readRDS(file.path(base_dir, 'transcript_pipeline_results.rds'))
gene_pipeline_results <-readRDS(file.path(base_dir, 'gene_pipeline_results.rds'))
tcc_pipeline_results <- readRDS(file.path(base_dir, 'tcc_pipeline_results.rds'))
tcc_table <- readRDS(file.path(base_dir, 'tcc_table.rds'))

print('reformatting tables')
gene_pipeline_results <- mutate(gene_pipeline_results, genes = target_id, gene_pval = pval, gene_qval = qval)
gene_pipeline_results$genes <- gsub('\\..*', '', gene_pipeline_results$genes)
gene_pipeline_results <- select(gene_pipeline_results, genes, gene_pval, gene_qval)


#convert ENSMUG names with period to that without
results$genes <- gsub('\\..*', '', results$genes)
results <- mutate(results, tx_lan_pval = lan, tx_min_pval = min, tx_lan_qval = lan.adjust, tx_min_qval = min.adjust)
results <- select(results, genes, tx_lan_pval, tx_min_pval, tx_lan_qval, tx_min_qval)

tcc_pipeline_results <- tcc_pipeline_results[order(tcc_pipeline_results$genes),]
tcc_pipeline_results <- mutate(tcc_pipeline_results, tcc_pval = pval, tcc_qval = qval)
tcc_pipeline_results <- select(tcc_pipeline_results, genes, tcc_pval, tcc_qval)

print('final_merge')
results <- merge(results, gene_pipeline_results, by = 'genes', all.x = TRUE)
results <- merge(results, tcc_pipeline_results, by ='genes', all.x = TRUE)
results <- results[ order(results$genes), ]

print('calculating weights')
weights_table <- tcc_table %>% group_by(genes) %>% summarise(weight = sum(meancounts))
weights_table <- weights_table[order(weights_table$genes),]
weights_table <- weights_table[-1,]
weights_table$genes <- as.character(weights_table$genes)

if(!all.equal(weights_table$genes, results$genes))
{
	print('ERROR: not all equal')
}

get_tgd_genes <- function(gene_names, qval)
{
	table <- data.frame(genes = gene_names, qval = qval)
	table$qval[is.na(table$qval)] <- 1
	#choose not to remove the filtered genes in order to use full transcriptome
	#table <- filter(table, !is.na(qval))
	genes <- as.factor(as.numeric(table$qval < .05))
	names(genes) <- as.character(table$genes)
	genes
}


select_top <- function(gene_names, qval, num)
{
	table <- data.frame(genes = gene_names, qval = qval)
	table$qval[is.na(table$qval)] <- 1
	table <- table[order(table$qval),]
	selected_genes <- table$genes[1:num]
	genes <- as.factor(as.numeric(table$genes %in% selected_genes))	
	names(genes) <- as.character(table$genes)
	genes
}

do_classical_analysis <- function(gene_names, qval, method_name=NULL)
{
	#nodeSize = minimum size of node as cutoff for consideration
	tgd <- new( "topGOdata", ontology='BP', allGenes = get_tgd_genes(gene_names, qval), nodeSize=5,
		   annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl" )
	resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
	tgd_results <- GenTable( tgd, Fisher.classic = resultTopGO.classic, orderBy = "Fisher.classic" , topNodes = 8390)
	tgd_results$Fisher.classic <- as.numeric(tgd_results$Fisher.classic)
	tgd_results$Fisher.classic[is.na(tgd_results$Fisher.classic)] <- 10^-30
	tgd_results$Fisher.classic[which(tgd_results$Significant < tgd_results$Expected)] <- 1
	tgd_results$method <- method_name
	tgd_results <- mutate(tgd_results, fwer = Fisher.classic)
	tgd_results
}


#agg analysis against classical
do_agg_analysis <- function(gene_names, pval, weight, go_map = NULL, method_name = NULL) 
{
	results <- data.frame(genes = gene_names, pval = pval, weight = weight)
	all_genes <- get_tgd_genes(unique(gene_names), c(0, rep(1, length(unique(gene_names))-1)))
	tgd <- new( "topGOdata", ontology='BP', allGenes = all_genes, nodeSize=5,
	   annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl" )
	git <- genesInTerm(tgd)
	print('performing aggregation')
	my_go_pval <- sapply(seq_along(git), function(i) map_GO(results, git[[i]]))
	fwer <- p.adjust(my_go_pval, method='bonferroni')
	tgd_results <- data.frame(GO.ID = names(git), my_go_pval = my_go_pval, fwer = fwer)
	tgd_results$method <- method_name
	tgd_results <- merge(tgd_results, go_map, by='GO.ID', all.x = TRUE)
	tgd_results
}

#this table gives annotation names
map_GO<- function(results_table, go_genes)
{
	table <- filter(results_table, genes %in% go_genes)
	lancaster(table$pval, table$weight)
}

print('doing classical analyses')
lan_tx <- do_classical_analysis(results$genes, results$tx_lan_qval, 'Lancaster Tx')
min_tx <- do_classical_analysis(results$genes, results$tx_min_qval, 'Sidak Tx')
gene <- do_classical_analysis(results$genes, results$gene_qval, 'Gene')
tcc <- do_classical_analysis(results$genes, results$tcc_qval, 'TCC')

go_map <- select(lan_tx, GO.ID, Term)
print('doing new analyses')
lan_tx_new <- do_agg_analysis(results$genes, results$tx_lan_pval, weights_table$weight, go_map, 'Aggregate GO Lancaster Tx')
tcc_new <- do_agg_analysis(results$genes, results$tcc_pval, weights_table$weight, go_map, 'Aggregate GO Lancaster TCC')
#tcc_direct <- do_agg_analysis(tcc_table$genes, tcc_table$sleuth_pval, weight = tcc_table$meancounts)
gene_new <- do_agg_analysis(results$genes, results$gene_pval, weights_table$weight, go_map, 'Aggregate GO Gene')
min_tx_new <- do_agg_analysis(results$genes, results$tx_min_pval, weights_table$weight, go_map, 'Aggregate GO Min Tx')

#plot -logpvalue of GO terms, filtered by name and minimum pvalue
plot_go_results <- function(go_results, term)
{
	go_results <- lapply(go_results, function(x)
		{
			print(x$method[1])
			print(colnames(x))
			y <- x[grep(term, x$Term), ]
			y <- dplyr::select(y, GO.ID, Term, fwer, method)
			y$logp <- -log(y$fwer, base=10)
			y
		})
	go_results <- do.call(rbind, go_results)
	
	filter_func <- function(go_results)
	{
		go_ids <- unique(go_results$GO.ID)
		in_vector <- sapply(go_ids, function(go_id)
			   {
					x <- filter(go_results, GO.ID == go_id)
			   		if(any(x$fwer < .05))
			   		{
			   				return(go_id)
			   		}
			   		NA
			   	})
		in_vector
	}
	
	filter_ids <- filter_func(go_results)
	filter_ids <- filter_ids[!is.na(filter_ids)]
	go_results <- filter(go_results, GO.ID %in% filter_ids)
	
	p <- ggplot(go_results, aes(Term, logp)) + geom_bar(aes(fill=method), position='dodge', stat='identity')
	p <- p + theme(axis.text.x=element_text(angle=-90, hjust=0))
	p <- p + geom_hline(yintercept = -log(.05, base=10))
	p <- p + annotate('text', go_results$Term[[5]], -log(.1, base=10), label='0.05 FWER')
	p <- p + xlab('GO Term')
	p <- p + ylab('-log(p-value)')
	png('~/GO.png', width = 10000, height = 7000, res=1200)
	print(p)
	dev.off()
}

#fisher exact test to test enrichment of GO terms in GO analysis
calculate_enrichment <- function(results, term)
{
	terms <- grep(term,results$Term)
	terms <- results[terms,]
	n_found <- sum(results$fwer < .05, na.rm=TRUE)
	n_total <- nrow(results)
	n_found_terms <-  sum(terms$fwer <.05, na.rm=TRUE)
	n_total_terms <- nrow(terms)
	a = n_found_terms
	b = n_found - n_found_terms
	c = n_total_terms - n_found_terms
	d = n_total - a - b - c
	x <- matrix(c(a,b,c,d), 2)
	p <- fisher.test(x, alternative = "greater")
	print(x)
	p
}

plot_scatter <- function(table1, table2)
{
	data <- merge(table1, table2, by = 'GO.ID')
	print(all.equal(data$Term.x, data$Term.y))
	x <- grep('immune', data$Term.x)
	terms <- as.numeric(1:nrow(data) %in% x)
	data$Terms <- as.factor(terms)
	data <- data[order(data$Terms),]
	png('~/scattergo.png', width=10000, height = 8000, res = 1200)
	p <- ggplot() + geom_point(data = data, aes(x=fwer.x, y = fwer.y, colour=Terms)) + 
			scale_colour_manual(values=c(rgb(0,0,0,.1), rgb(1,0,0,.7)),
								labels = c('GO Terms without "immune"', 
										   'GO Terms with "immune"'))
	p <- p + xlab(paste0('Enrichment Test (p-value)')) + ylab(paste0('Perturbation Test (p-value)'))
	print(p)
	dev.off()
}


print('plotting go results')
go_results <- list(lan_tx, min_tx, gene, tcc)
plot_go_results(go_results, 'immune')

print('plot scatter')


go_results <- list(tcc, tcc_new)
enrichments <- lapply(go_results, function(x) calculate_enrichment(x, 'immune'))
enp <- sapply(enrichments_1, function(x) x$p.value)

