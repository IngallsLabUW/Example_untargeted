# Could Cloud Computing?
#   wkumler
# August 25, 2021
# Cloud computing for XCMS parameter optimization
# Goal: compare lots of different XCMS peak-picking parameters to the “gold standard” of manually-integrated peaks
# Data: 24 Falkor samples, integrated by Raafay?
#   Params to alter:
#   The main parameters that I’m interested in are all in the centWaveParam object, from and for XCMS:

library(xcms)


CentWaveParam()

# I’ve fiddled with a lot of the above settings but never systematically, 
# and they can all have pretty major effects on peakpicking quality.
# What I imagine using cloud computing for is iterating over lots of different parameters 
# to see how well they compare to the gold standard.
# This would essentially be done by running XCMS a bunch of times with different parameters, 
# then annotating the peaklist and calculating some kind of correlation between 
# the calculated areas and the manual integrations.

library(tidyverse)
param_df <- expand_grid(
  ppm=c(1, 2.5, 5, 10, 20),
  peakwidth_min=c(5, 10, 15, 25, 50),
  peakwidth_max=c(5, 10, 15, 25, 50),
  integrate=c(1, 2),
  qscore=c(5, 10, 20, 50),
  extendLengthMSW=c(TRUE, FALSE)
) %>%
  filter(peakwidth_min<=peakwidth_max)

# Then, run XCMS in some kind of loop, providing different parameters each time:
ms_files <- list.files("~/work/untargeted/mzMLs/neg", pattern = "mzML")
raw_data <- readMSData(ms_files, mode = "onDisk")
ingalls_stans <- read.csv("http://github/path/to/regina/standards")
raafay_vals <- read.csv("path/to/raafay/integrations.csv")

for(row in seq_len(nrow(param_df))){
  cwp <- CentWaveParam(
    ppm=param_df$ppm[row],
    integrate=param_df$integrate[row],
    qscore=param_df$qscore[row],
    peakwidth=c(param_df$peakwidth_min[row],
                param_df$peakwidth_max[row])
  )
  xdata <- findChromPeaks(raw_data, param = cwp)
  #...other XCMS steps & Will custom code
  xdata_peaks <- chromPeakData(xdata) %>%
    filter(feature_name%in%ingalls_stans$compound_name)
  similarity_vals <- sapply(unique(xdata_peaks$feature_name), function(j){
    xcms_integrations <- xdata_peaks[xdata_peaks$feature_name==j,]
    raafay_integrations <- raafay_vals[raafay_vals$feature_name==j,]
    cor(xcms_integrations, raafay_integrations)
  })
}
# The question is how to best parallelize the process, 
# which is dependent upon the cloud’s architecture. 
# We can parallelize by file, which is how we’ve been doing it on laptops and the PLGS:
register(BPPARAM = SnowParam(tasks = length(ms_files)))
# All the above things or we can do a single file per core, 
# which would probably make best use of our computing hours 
# (assuming each set of params takes about the same amount of time to compute) 
# but I have no idea how this should be done and will depend on resources 
# available to each node (memory, storage, floprate)

#!/bin/bash
# Job name:
#SBATCH --job-name=XCMSParams_detailed_params_here
#
# Account:
#SBATCH --account=whatever
#
# Partition:
#SBATCH --partition=something_from_cloud
#
# Request one node:
#SBATCH --nodes=number_of_param_permutations_to_run
#
# Specify number of tasks for use case:
#SBATCH --ntasks-per-node=number_of_files (maybe?)
#
# Processors per task:
#SBATCH --cpus-per-task=number_of_files (maybe?)
#
# Wall clock limit:
#SBATCH --time=02:00:00
#
# Log file location:
#SBATCH --output="../SLURMlogs/job_%j.out"


## Command(s) to run:
# module load r/3.5.1
# module load r-packages/default
# module load r-packages/custom_ingalls
# 
# R CMD BATCH rscripts/XCMS_params_given_params.R