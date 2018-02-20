filter_func <- function(row)
{
	mean(row) > 0
}

normalize_tccs <- function(counts)
{
		counts <- log(counts + 1)
		counts <- counts / sum(counts)
		counts
}

DESeq_wt <- function(counts, design)
{
	#need nsamples =  ncol(count) = nrow(design)
	DESeqData <- DESeqDataSetFromMatrix(counts, DataFrame(design), ~ factor(design))
	estimateSizeFactors(DESeqData) -> DESeqData
	estimateDispersions(DESeqData) -> DESeqData
	nbinomWaldTest(DESeqData, betaPrior=FALSE)-> DESeqData
	DESeqData
}

DESeq_lrt <- function(counts, design)
{
	DESeqData <- DESeqDataSetFromMatrix(counts, DataFrame(design), ~ factor(design))
	DESeqData <- estimateSizeFactors(DESeqData)
	DESeqData <- DESeq2::DESeq(DESeqData, test = 'LRT', reduced = ~1)
	DESeqData
}

sleuth_wt_pipeline <- function(so)
{
	so <- sleuth_fit(so, ~condition, 'full')
	so <- sleuth_wt(so, 'condition1', 'full')
	sleuth_table <- sleuth_results(so, 'condition1', which_model = 'full')
	sleuth_table
}

sleuth_lrt_pipeline <- function(so)
{
	so <- sleuth_fit(so, ~condition, 'full')
	so <- sleuth_fit(so, ~1, 'null')
	so <- sleuth_lrt(so, 'null', 'full')
	sleuth_table <- sleuth_results(so, 'null:full', 'lrt')
	sleuth_table
}

limma <- function(counts, design)
{
	design <- model.matrix(~design)
	v <- voom(counts, design)
	vfit <- lmFit(v, design)
	efit <- eBayes(vfit)
	#tfit <- treat(efit, lfc=1)
	#tfit
	efit
}

#map = map from tcc to transcript id
#transcripts = map from transcript to genes
make_tcc_to_gene_map <- function(transcripts, map)
{
	ntccs <- nrow(map)
	tcc2gene <- lapply(1:ntccs, function(i) map_to_gene(i, transcripts, map))
	saveRDS(tcc2gene, './tcc2gene.rds')
	tcc2gene
}


make_tcc_to_gene_map2 <- function(transcripts, map)
{
	ntccs <- nrow(map)
	tcc2gene <- sapply(1:ntccs, function(i) map_to_gene2(i, transcripts, map))
	tcc2gene
}

#hidden, called by make_tcc_to_gene_map
map_to_gene2 <- function(i, transcripts, map)
{
	#map tcc i to transcript id
	transcript_indices <- unlist(strsplit(toString(map[[i,2]]), ','))
	transcript_indices <- sapply(transcript_indices, strtoi)
	#have to map from 0 indexed transcript indices to the t2g map
	transcript_indices <- transcript_indices +1
	g <- sapply(transcript_indices, function(x) transcripts[x,2])
	unique_g <- unique(g)
	
	if(length(unique_g) == 1)
	{
		return(unique_g)
	}
	else
	{
		return('')
	}
}


#hidden, called by make_tcc_to_gene_map
map_to_gene <- function(i, transcripts, map)
{
	#map tcc i to transcript id
	transcript_indices <- unlist(strsplit(toString(map[[i,2]]), ','))
	transcript_indices <- sapply(transcript_indices, strtoi)
	#have to map from 0 indexed transcript indices to the t2g map
	transcript_indices <- transcript_indices +1
	t <- sapply(transcript_indices, function(x) transcripts[x,1])
	g <- sapply(transcript_indices, function(x) transcripts[x,2])
	unique_g <- unique(g)
	list(tcc_id=i, t=t, g=g, unique = unique_g)
}

#TODO factor this out
make_tcc_table <- function(tcc2genemap)
{
	n <- length(tcc2genemap)
	#df <- data.frame(pvalues=rep(1.0, n), genes = character(n), stringsAsFactors=FALSE)
	#for (i in 1:n)
	#{
	#	if(length(tcc2gene[[i]]$unique) == 1)
	#	{
	#		df[i,2] <- tcc2gene[[i]]$unique
	#	}
	#}
	#df

	table <- sapply(1:n, function(i) check_unambiguous(tcc2genemap[[i]]))
	table <- unlist(table)
	table <- matrix(table, ncol=2, byrow=TRUE)
	table <- as.data.frame(table)
	table[,1] <- as.numeric(table[,1])
	colnames(table) <- c('pvalue', 'genes')
	table
}

check_unambiguous <- function(tcc)
{
	if(length(tcc$unique)==1)
	{
		return(list(1.0, tcc$unique))
	}
	return(list(1.0, ''))
}

filter_all_nonunique <- function(tcc2gene, genes)
{
	genes <- unique(genes)
	genes <- genes[order(genes)]
	genes_table <- data.frame(genes = genes, unique = TRUE, stringsAsFactors=FALSE)
	print(head(genes_table))
	n <- length(tcc2gene)
	for (i in 1:n)
	{
		if(length(tcc2gene[[i]]$unique) > 1)
		{
			for(j in 1:length(tcc2gene[[i]]$unique))
			{
				gene <- tcc2gene[[i]]$unique[j]
				genes_table$unique[which(genes_table$genes == gene)] <- FALSE
			}

		}
	}
	genes_table
}

