
MinMethod <- function(pvalues)
{
	if(length(pvalues) == 0)
	{
			return(NA)
	}
	pvalues <- pvalues[!is.na(pvalues)]
	if(length(pvalues) == 0)
	{
			return(NA)
	}
	n <- length(pvalues)
	m <- min(pvalues)
	1 - (1-m) ^ n
}

DESeq <- function(counts, design)
{
	#need nsamples =  ncol(count) = nrow(design)
	DESeqData <- DESeqDataSetFromMatrix(counts, DataFrame(design), ~ factor(design))
	#DESeqData <- DESeq(DESeqData, test=c("Wald"), fitType='parametric', quiet=FALSE)
	estimateSizeFactors(DESeqData) -> DESeqData
	estimateDispersions(DESeqData) -> DESeqData
	nbinomWaldTest(DESeqData, betaPrior=FALSE)-> DESeqData

	DESeqData
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

make_tcc_table <- function(tcc2genemap)
{
	n <- length(tcc2genemap)
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


fisher <- function(pvalues)
{
	pvalues <- pvalues[!is.na(pvalues)]
	if(length(pvalues)==0)
	{
		return (NA)
	}
	chisq = -2 * sum(log(pvalues))
	df <- 2* length(pvalues)
	pchisq(chisq, df, lower.tail = FALSE)
}

lancaster <- function(pvalues, weights)
{
	weights <- weights[!is.na(pvalues)]
	pvalues <- pvalues[!is.na(pvalues)]
	pvalues <- pvalues[weights>0]
	weights <- weights[weights>0]
	
	if(length(weights) != length(pvalues))
	{
		print('error, weights not equal to pvalues')
	}

	if(length(pvalues) == 0)
	{
		return(NA)
	}
	if(!any(weights))
	{
		return(NA)
	}
	if(length(pvalues) == 1)
	{
		return(pvalues)
	}
	t <- sapply(1:length(pvalues), function(i) lts(pvalues[i], weights[i]))
	t <- sum(t)
	p <- pchisq(t, sum(weights), lower.tail=FALSE) 
	return(p)
}

lts <- function(pvalue, weight)
{
	qgamma(pvalue, shape = weight /2, scale = 2, lower.tail=FALSE)
}



