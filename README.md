This repository contains the scripts and software for reproducing the results and figures of the paper "Gene-level differential analysis at transcript-level resolution" by Lynn Yi, Harold Pimentel, Nicolas L Bray and Lior Pachter. The code can also be used to apply the aggregation methods described in the paper to new datasets. The software in the repository was written by Lynn Yi. 

R/Snakefile is an example pipeline for downloading fastq files, performing pseudoquant, and bootstraping on TCCs. The remainder processes for calling sleuth and aggregation p-values are performed in R scripts, tcc_pipeline.R and transcript_pipeline.R

R/aggregation.R contains logic for performing aggregation, incuding mapping TCCs to genes.
R/tcc2bootstrap.R contains logic for performing bootstraps on TCCs and writing h5 files that sleuth can take as input.

The folders SRPXXXXX contain code for reproducing analysis for the two datasets in the paper. They include Snakefiles for read downloading and quantification,  aggregation pipelines, and GO analysis. plot_transcripts.R include code for reproducing Figures 1 and 2 in paper. topGO.R include code for reproducing Figure 5 and performing GO analysis.

The folder simulation_pipeline contains code for reproducing the analyses of simulations described in the paper. pachterlab/sleuthpaperanalysis must be utilized first to create simulations. Then simulation_pipeline.R will run various differential expression and aggregation methods. averaging_fdrs.R and roc_curve.R will handle averaging FDR and sensitivities. Finally mamabear is invoked in mamabear.R to plot.

