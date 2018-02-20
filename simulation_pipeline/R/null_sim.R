#simulated p-values under the null

library(dplyr)
calculate_fdr_min <- function(pval, pvalues, n)
{
	if(length(pvalues) != length(n))
	{
		print('there is ERROR in length')
    }
    if(is.na(pval))
    {
   		return(NA)
    }
   numerator <- sum(1 - (1-pval)^n)
   denom <- sum(pvalues <= pval)
   fdr <- numerator / denom
   fdr
}


t2g <- read.table('~/transcriptomes/Mus_musculus.GRCm38.cdna.rel.88.transcripts', sep='\t')
table <- t2g %>% group_by(V2) %>% summarise(n=n())
n <- table$n
pvals <- lapply(n, function(n) runif(n, min = 0, max=1))
min_p <- sapply(pvals, min)
sidak <- sapply(seq_along(pvals), function(i) {1 -  (1 - min(pvals[[i]]))^n[i]})

results <- data.frame(n, min_p, sidak)
p <- sapply(results$min_p, function(x) calculate_fdr_min(x, results$min_p, results$n))

