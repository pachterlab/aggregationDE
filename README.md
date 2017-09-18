Written by Lynn Yi on Sep 17, 2017

This repository contains the scripts and software for reproducing the results and figures of the paper "Gene-level differential analysis at transcript-level resolution" by Lynn Yi, Harold Pimentel, Nicolas L Bray and Lior Pachter. The code can also be used to apply the aggregation methods described in the paper to new datasets. The software in the repository was written by Lynn Yi. 

R/Snakefile is an example pipeline for downloading fastq files, performing pseudoalignment, and bootstraping on TCCs. The remaining processes for calling sleuth and aggregating p-values are performed in the R scripts tcc_pipeline.R and transcript_pipeline.R

R/aggregation.R contains code for performing aggregation, incuding mapping TCCs to genes. R/tcc2bootstrap.R contains code for performing bootstraps on TCCs and writing h5 files that sleuth can take as input.

The folders SRPXXXX contain code for reproducing analysis for the two biological datasets analyzed in the paper. They include Snakefiles for downloading reads from the short read archive and for quantification, aggregation pipelines, and GO analysis.

The folder simulation_pipeline contains code for reproducing the analyses of simulations described in the paper.