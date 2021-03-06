
#pipeline to "bootstrap" from TCCs
#output into kallisto abundance.h5 format in order for sleuth to read and analyze it

#assume df is TCC matrix
#assume h5template is template

library(rhdf5)

args <- commandArgs(TRUE)
experiment <- args[[1]]
exp_string <- args[[2]]
print(experiment)
print(exp_string)

tcc_path <- file.path(experiment, exp_string, 'matrix.tsv')
directory <- file.path(experiment, exp_string, 'sleuth_tcc')

tccs <- scan(tcc_path, skip=1, what = numeric(), sep='\t')
tccs <- matrix(tccs, nrow=7)
tccs <- tccs[-1,]
print(dim(tccs))

ncells <- nrow(tccs)
ntccs <- ncol(tccs)
print('ncells, ntccs')
print(ncells)
print(ntccs)

write_h5 <- function(counts, h5file)
{
	print('creating file')
	print(h5file)
	h5createFile(h5file)
	h5createGroup(h5file, 'aux')
	h5createGroup(h5file, 'bootstrap')
	h5write(counts, h5file, 'est_counts')
	h5write(1:ntccs, h5file, 'aux/ids')
	h5write(rep(1.0, ntccs), h5file, 'aux/eff_lengths') 
	h5write(rep(1.0, ntccs), h5file, 'aux/lengths') 
	h5write(rep(0, 1000), h5file, 'aux/fld') 
	h5write(rep(1, 4096), h5file, 'aux/bias_observed')
	h5write(rep(1, 4096), h5file, 'aux/bias_normalized')
	h5write(30, h5file, 'aux/num_bootstrap')
	h5write(10, h5file, 'aux/index_version')
	h5write(42.5, h5file, 'aux/kallisto_version')
	h5write(as.character(Sys.time()), h5file, 'aux/start_time')
	print('writing bootstraps')	
	bootstraps <- rmultinom(30, sum(counts), counts)
	sapply(1:30, function(i) h5write(bootstraps[,i], h5file, paste0('bootstrap/bs',(i-1))))
	H5close()
}

sapply(1:ncells, function(i) write_h5(tccs[i,], file.path(directory, i, 'abundance.h5')))
print('finished making bootstraps for tccs')

