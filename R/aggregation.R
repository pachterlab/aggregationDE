

make_tcc_to_gene_map <- function(transcripts, map)
{
	ntccs <- nrow(map)
	tcc2gene <- sapply(1:ntccs, function(i) map_to_gene(i, transcripts, map))
	tcc2gene
}

map_to_gene <- function(i, transcripts, map)
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



