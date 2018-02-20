
library(tximport)

x <- list.files('/home/lynnyi/tximport/R/')
y <- file.path('/home/lynnyi/tximport/R/', x)
sapply(y, source)
print('sourcing')
print(y)

base_path <- '/home/lynnyi/SRP100701'
#s2c <- read.table(file.path(base_path,'sample_table.txt'), header=TRUE, sep='\t')

#samples <- s2c$Run_s
samples <- c(condition1, condition2)

files <- sapply(samples, function(x) file.path(base_path, 'kallisto', x, 'abundance.h5'))
transcripts <- read.table('~/transcriptomes/Mus_musculus.GRCm38.cdna.rel.88.transcripts')
tximp <- tximport(files, type = 'kallisto', tx2gene = transcripts)
gene_counts <- round(tximp$counts)



