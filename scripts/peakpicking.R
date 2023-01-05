# Script to run XCMS peakpicking, retention time correction, and grouping
# Called by Control.Rmd
# Requires dev version of XCMS for improved peakpicking given settings

# ms_files comes from Control.Rmd
# qscore_threshold comes from Control.Rmd


# Read in the raw data ----
message("Reading files...")
start_time <- Sys.time()
raw_data <- readMSData(files = normalizePath(paste0(mzml_path, ms_files)), 
                       msLevel. = 1, centroided. = TRUE, mode = "onDisk",
                       pdata = as(cruise_metadata, "AnnotatedDataFrame"))
message("Time to read in files: ", round(Sys.time()-start_time, digits = 2), " min")

# Perform peakpicking ----
message("Picking peaks...")
start_time <- Sys.time()
xdata <- suppressMessages(findChromPeaks(raw_data, param = cwp))
message("Time to perform peakpicking: ", 
        round(Sys.time()-start_time, digits = 2), " min")
# 13 minutes


# Assign new quality scores ----
# speedyQscoreCalculator is custom function which should be in functions script
# sourced by Control.Rmd
message("Assigning quality scores...")
start_time <- Sys.time()
xcms_peakdf <- chromPeaks(xdata) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  mutate(filename=fileNames(xdata)[.[["sample"]]]) %>%
  arrange(mz)
split_xcms_filedata <- split(xcms_peakdf, xcms_peakdf$filename)
files_qscores <- bplapply(X = split_xcms_filedata, FUN = speedyQscoreCalculator,
                          grabSingleFileData, qscoreCalculator)

peakdf_qscored <- files_qscores %>%
  do.call(what = rbind) %>%
  as.data.frame(stringsAsFactors=FALSE) %>% 
  filter(sn>qscore_threshold) %>%
  select(-filename) %>%
  arrange(sample, rtmin, rtmax) %>%
  as.matrix()
xdata_cleanpeak <- `chromPeaks<-`(xdata, peakdf_qscored)
message("Time to assign quality scores: ", 
        round(Sys.time()-start_time, digits = 2), " min")
# 20 minutes



# Other XCMS things (rtcor, group) ----
message("Other XCMS things:")
start_time <- Sys.time()
message("Retention time correction (double progress bar)...")
xdata_rt <- suppressMessages(adjustRtime(xdata_cleanpeak, param = obp))
# plotAdjustedRtime(xdata_rt)

message("Grouping...")
xdata_cor <- groupChromPeaks(xdata_rt, param = pdp)

message("Filling peaks...")
xdata_filled <- suppressMessages(fillChromPeaks_wkumler(xdata_cor, param = fpp))

feature_defs <- featureDefinitions(xdata_filled)
max_features <- unique(nchar(rownames(feature_defs)))-2
raw_peaks <- lapply(seq_len(nrow(feature_defs)), function(i){
  cbind(feature=sprintf(paste0("FT%0", max_features, "d"), i),
        peak_id=unlist(feature_defs$peakidx[i]))
}) %>%
  do.call(what=rbind) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  mutate(peak_id=as.numeric(peak_id)) %>%
  cbind(chromPeaks(xdata_filled)[.$peak_id, ]) %>%
  mutate(filename=basename(fileNames(xdata_filled))[sample]) %>%
  arrange(feature, sample) %>%
  mutate(across(starts_with("rt"), round, digits=5))

message("Time for other XCMS things: ", 
        round(Sys.time()-start_time, digits = 2), " min")
# 30 minutes
