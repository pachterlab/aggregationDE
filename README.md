
Written by Lynn Yi on Sep 17, 2017

AggregationDE is a pipeline for aggregation p-values of transcripts and TCCs to obtain gene-level differential expression.

R/Snakefile is an example pipeline for downloading fastq files, performing pseudoquant, and bootstraping on TCCs. The remainder processes for calling sleuth and aggregation p-values are performed in R scripts, tcc_pipeline.R and transcript_pipeline.R

R/aggregation.R contains logic for performing aggregation, incuding mapping TCCs to genes.
R/tcc2bootstrap.R contains logic for performing bootstraps on TCCs and writing h5 files that sleuth can take as input.

The folders SRPXXXX contain code for reproducing analysis for the two datasets in the paper. They include Snakefiles for read downloading and quantification,  aggregation pipelines, and GO analysis.

The folder simulation_pipeline contains code for reproducing analysis from simulations.

