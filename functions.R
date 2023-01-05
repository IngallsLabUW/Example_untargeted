# Functions called by Control.Rmd or scripts within
# Should be sourced every time!

cv <- function(x)sd(x)/mean(x)

grabSingleFileData <- function(filename){
  msdata <- mzR:::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  spectra_list <- lapply(seq_len(nrow(fullhd)), function(x){
    given_peaks <- mzR::peaks(msdata, x)
    rtime <- fullhd[x, "retentionTime"]
    return(cbind(rtime, given_peaks))
  })
  all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)), 
                        c("rt", "mz", "int"))
  return(all_data)
}
speedyQscoreCalculator <- function(file_peaks, grabSingleFileData, qscoreCalculator){
  library(data.table)
  file_data <- grabSingleFileData(unique(file_peaks$filename))
  file_data_table <- as.data.table(file_data)
  peak_splits <- split(file_peaks, ceiling(seq_len(nrow(file_peaks))/100))
  file_qscores <- lapply(peak_splits, function(i){
    eic_many <- file_data_table[mz%between%c(min(i$mzmin), max(i$mzmax))]
    individual_peaks <- split(i, seq_len(nrow(i)))
    output <- sapply(individual_peaks, function(peak_row){
      eic <- eic_many[rt%between%c(peak_row$rtmin, peak_row$rtmax) & 
                        mz%between%c(peak_row$mzmin, peak_row$mzmax)]
      qscoreCalculator(eic)
    }, USE.NAMES = FALSE)
  })
  file_peaks$sn <- unlist(file_qscores)
  return(file_peaks)
}
qscoreCalculator <- function(eic){
  #Check for bogus EICs
  if(nrow(eic)<5){
    return(0)
  }
  #Calculate where each rt would fall on a beta dist (accounts for missed scans)
  scaled_rts <- (eic$rt-min(eic$rt))/(max(eic$rt)-min(eic$rt))
  
  # Create a couple different skews and test fit
  maybe_skews <- c(2.5,3,4,5) #Add 7 to catch more multipeaks and more noise
  #Add 2 to catch very slopey peaks and more noise
  best_skew <- maybe_skews[which.max(sapply(maybe_skews, function(x){
    cor(dbeta(scaled_rts, shape1 = x, shape2 = 5), eic$int)
  }))]
  perf_peak <- dbeta(scaled_rts, shape1 = best_skew, shape2 = 5)
  peak_cor <- cor(perf_peak, eic$int)
  
  
  #Calculate the normalized residuals
  residuals <- eic$int/max(eic$int)-perf_peak/max(perf_peak)
  #Calculate the minimum SD, after normalizing for any shape discrepancy
  old_res_sd <- sd(residuals)
  norm_residuals <- diff(residuals)
  new_res_sd <- sd(norm_residuals)
  while(new_res_sd<old_res_sd){
    old_res_sd <- new_res_sd
    norm_residuals <- diff(residuals)
    new_res_sd <- sd(residuals)
  }
  #Calculate SNR
  SNR <- (max(eic$int)-min(eic$int))/sd(norm_residuals*max(eic$int))
  #Return the quality score
  return(round(SNR*peak_cor^4*log10(max(eic$int))))
}
qscoreCalculator <- compiler::cmpfun(qscoreCalculator)
isIsoAdduct <- function(file_peaks, xdata, grabSingleFileData, 
                        checkPeakCor, pmppm, trapz, polarity, mzml_path){
  file_peaks <- file_peaks[order(file_peaks$rtmax),]
  file_path <- paste0(mzml_path, unique(file_peaks$filename))
  file_data <- grabSingleFileData(file_path)
  file_data$rt <- xcms::adjustedRtime(xdata)[
    MSnbase::fromFile(xdata)==unique(file_peaks$sample)][factor(file_data$rt)]
  file_data$rt <- round(file_data$rt, digits = 5)
  library(data.table)
  file_data_table <- as.data.table(file_data)
  
  peak_splits <- split(file_peaks, ceiling(seq_len(nrow(file_peaks))/10))
  
  iso_matches_all <- lapply(peak_splits, function(i){
    eic_many <- file_data_table[rt>min(i$rtmin)&rt<max(i$rtmax)]
    individual_peaks <- split(i, seq_len(nrow(i)))
    iso_matches <- lapply(individual_peaks, function(peak_row_data){
      isos_to_check <- c(C13=-1.003355, X2C13=-2*1.003355, S34=-1.995796, S33=-0.999387,
                         N15=-0.997035, O18=-2.004244) + peak_row_data[["mz"]]
      if(polarity=="pos"){
        adducts_to_check <- c(Na=-22.98922+1.007276, NH4=-18.0338+1.007276,
                              mH2OpH=+18.0106, K=-38.963708+1.007276) + peak_row_data[["mz"]]
        more_adducts <- c(X2H=peak_row_data[["mz"]]*2-2*1.007276+1.007276)
      } else if(polarity=="neg"){
        adducts_to_check <- c(Cl=-34.969402+1.007276, Ac=-59.013851+1.007276,
                              mH2OmH=+18.0106, Br=-78.918885+1.007276) + peak_row_data[["mz"]]
        more_adducts <- c(X2H=peak_row_data[["mz"]]*2-2*1.007276-1.007276)
      } else {
        stop("Unrecognized polarity in adduct search, see function isIsoAdduct.")
      }
      masses_to_check <- c(isos_to_check, adducts_to_check, more_adducts)
      
      init_eic <- eic_many[mz%between%pmppm(peak_row_data["mz"], ppm = 5) & 
                             rt%between%c(peak_row_data["rtmin"], peak_row_data["rtmax"])]
      
      output <- lapply(masses_to_check, checkPeakCor, rtmin=peak_row_data["rtmin"], 
                       rtmax=peak_row_data["rtmax"], init_eic = init_eic, 
                       file_data_table = eic_many, pmppm = pmppm, trapz = trapz)
      linear_output <- do.call(c, output)
      names(linear_output) <- paste0(rep(names(masses_to_check), each=2), 
                                     c("_match", "_area"))
      if(!nrow(init_eic)){
        init_area <- 0
      } else {
        init_area <- trapz(init_eic$rt, init_eic$int)
      }
      linear_output <- c(M_area=init_area, linear_output)
      return(linear_output)
    })
    iso_matches <- do.call(rbind, iso_matches)
    return(iso_matches)
  })
  iso_matches_all <- do.call(rbind, iso_matches_all)
  return(cbind(file_peaks, as.data.frame(iso_matches_all)))
}
checkPeakCor <- function(mass, rtmin, rtmax, init_eic, 
                         file_data_table, pmppm, trapz){
  given_eic <- file_data_table[mz%between%pmppm(mass, ppm = 5) & 
                                 rt%between%c(rtmin, rtmax)]
  if(nrow(given_eic)<5){
    return(c(0, 0))
  }
  merged_eic <- merge(init_eic, given_eic, by="rt")
  if(nrow(merged_eic)<5){
    return(c(0, 0))
  }
  
  # plot(merged_eic$rt, scales::rescale(merged_eic$int.y))
  # lines(init_eic$rt, scales::rescale(merged_eic$int.x))
  
  peak_match <- cor(merged_eic$int.x, merged_eic$int.y)
  peak_area <- trapz(merged_eic$rt, merged_eic$int.y)
  return(c(peak_match, peak_area))
}
checkPeakCor <- compiler::cmpfun(checkPeakCor)
trapz <- function(x, y) {
  m <- length(x)
  xp <- c(x, x[m:1])
  yp <- c(numeric(m), y[m:1])
  n <- 2*m
  p1 <- sum(xp[1:(n-1)]*yp[2:n]) + xp[n]*yp[1]
  p2 <- sum(xp[2:n]*yp[1:(n-1)]) + xp[1]*yp[n]
  
  return(0.5*(p1-p2))
}
trapz <- compiler::cmpfun(trapz)
findIsoAdduct <- function(file_peaks, xdata, grabSingleFileData,
                          checkPeakCor, pmppm, trapz, polarity, mzml_path){
  file_peaks <- file_peaks[order(file_peaks$rtmax),]
  file_path <- paste0(mzml_path, unique(file_peaks$filename))
  file_data <- grabSingleFileData(file_path)
  file_data$rt <- xcms::adjustedRtime(xdata)[
    MSnbase::fromFile(xdata)==unique(file_peaks$sample)][factor(file_data$rt)]
  file_data$rt <- round(file_data$rt, digits = 5)
  library(data.table)
  file_data_table <- as.data.table(file_data)
  
  peak_splits <- split(file_peaks, ceiling(seq_len(nrow(file_peaks))/10))
  
  iso_matches_all <- lapply(peak_splits, function(i){
    eic_many <- file_data_table[rt>min(i$rtmin)&rt<max(i$rtmax)]
    individual_peaks <- split(i, seq_len(nrow(i)))
    iso_matches <- lapply(individual_peaks, function(peak_row_data){
      isos_to_check <- c(C13=1.003355, X2C13=2*1.003355, S34=1.995796, S33=0.999387,
                         N15=0.997035, O18=2.004244) + peak_row_data[["mz"]]
      if(polarity=="pos"){
        adducts_to_check <- c(Na=22.98922-1.007276, NH4=18.0338-1.007276,
                              H2O_H=-18.0106, K=38.963708-1.007276) + peak_row_data[["mz"]]
        more_adducts <- c(X2H=peak_row_data[["mz"]]-1.007276+2*1.007276)/2
      } else if(polarity=="neg"){
        adducts_to_check <- c(Cl=34.969402+1.007276, Ac=59.013851+1.007276,
                              H2O_H=18.0106, Br=78.918885+1.007276) + peak_row_data[["mz"]]
        more_adducts <- c(X2H=(peak_row_data[["mz"]]+1.007276)*2-2*1.007276)
      } else {
        stop("Unrecognized polarity in adduct search, see function findIsoAdduct.")
      }
      masses_to_check <- c(isos_to_check, adducts_to_check, more_adducts)
      
      init_eic <- eic_many[mz%between%pmppm(peak_row_data["mz"], ppm = 5) & 
                             rt%between%c(peak_row_data["rtmin"], peak_row_data["rtmax"])]
      
      output <- lapply(masses_to_check, checkPeakCor, rtmin=peak_row_data["rtmin"], 
                       rtmax=peak_row_data["rtmax"], init_eic = init_eic, 
                       file_data_table = eic_many, pmppm = pmppm, trapz = trapz)
      linear_output <- do.call(c, output)
      names(linear_output) <- paste0(rep(names(masses_to_check), each=2), 
                                     c("_match", "_area"))
      
      if(!nrow(init_eic)){
        init_area <- 0
      } else {
        init_area <- trapz(init_eic$rt, init_eic$int)
      }
      
      linear_output <- c(M_area=init_area, linear_output)
      return(linear_output)
    })
    iso_matches <- do.call(rbind, iso_matches)
    return(iso_matches)
  })
  iso_matches_all <- do.call(rbind, iso_matches_all)
  return(cbind(file_peaks, as.data.frame(iso_matches_all)))
}
fillChromPeaks_wkumler <- function(object, param){
  msLevel <- 1L
  BPPARAM <- bpparam()
  
  if (length(msLevel) != 1) 
    stop("Can only perform peak filling for one MS level at a time")
  if (!hasFeatures(object, msLevel = msLevel)) 
    stop("No feature definitions for MS level ", 
         msLevel, " present. Please run 'groupChromPeaks' first.")
  if (xcms:::.hasFilledPeaks(object))
    message("Filled peaks already present, adding still missing", 
            " peaks.")
  if (xcms:::hasChromPeaks(object) & !xcms:::.has_chrom_peak_data(object)) 
    object <- updateObject(object)
  startDate <- date()
  expandMz <- expandMz(param)
  expandRt <- expandRt(param)
  fixedMz <- fixedMz(param)
  fixedRt <- fixedRt(param)
  ppm <- ppm(param)
  message("Defining peak areas for filling-in .", 
          appendLF = FALSE)
  fdef <- featureDefinitions(object, msLevel = msLevel)
  aggFunLow <- median
  aggFunHigh <- median
  tmp_pks <- chromPeaks(object)[, c("rtmin", "rtmax", 
                                    "mzmin", "mzmax")]
  pkArea <- do.call(rbind, lapply(fdef$peakidx, function(z) {
    pa <- c(aggFunLow(tmp_pks[z, 1]), aggFunHigh(tmp_pks[z, 2]), 
            aggFunLow(tmp_pks[z, 3]), aggFunHigh(tmp_pks[z, 4]))
    if (ppm != 0) {
      mzmean <- mean(pa[3:4])
      tittle <- mzmean * (ppm/2)/1e+06
      if ((pa[4] - pa[3]) < (tittle * 2)) {
        pa[3] <- mzmean - tittle
        pa[4] <- mzmean + tittle
      }
    }
    if (expandRt != 0) {
      diffRt <- (pa[2] - pa[1]) * expandRt/2
      pa[1] <- pa[1] - diffRt
      pa[2] <- pa[2] + diffRt
    }
    if (expandMz != 0) {
      diffMz <- (pa[4] - pa[3]) * expandMz/2
      pa[3] <- pa[3] - diffMz
      pa[4] <- pa[4] + diffMz
    }
    if (fixedMz != 0) {
      pa[3] <- pa[3] - fixedMz
      pa[4] <- pa[4] + fixedMz
    }
    if (fixedRt != 0) {
      pa[1] <- pa[1] - fixedRt
      pa[2] <- pa[2] + fixedRt
    }
    pa
  }))
  rm(tmp_pks)
  message(".", appendLF = FALSE)
  colnames(pkArea) <- c("rtmin", "rtmax", "mzmin", "mzmax")
  pkArea <- cbind(group_idx = 1:nrow(pkArea), pkArea, mzmed = as.numeric(fdef$mzmed))
  pkGrpVal <- featureValues(object, value = "index", msLevel = msLevel)
  message(".", appendLF = FALSE)
  if (!any(is.na(rowSums(pkGrpVal)))) {
    message("No missing peaks present.")
    return(object)
  }
  message(".", appendLF = FALSE)
  objectL <- vector("list", length(fileNames(object)))
  pkAreaL <- objectL
  req_fcol <- requiredFvarLabels("OnDiskMSnExp")
  min_fdata <- fData(object)[, req_fcol]
  rt_range <- range(pkArea[, c("rtmin", "rtmax")])
  if (hasAdjustedRtime(object)) 
    min_fdata$retentionTime <- adjustedRtime(object)
  for (i in 1:length(fileNames(object))) {
    fd <- min_fdata[min_fdata$fileIdx == i, ]
    fd$fileIdx <- 1L
    objectL[[i]] <- new("OnDiskMSnExp", 
                        processingData = new("MSnProcess", 
                                             files = fileNames(object)[i]), 
                        featureData = new("AnnotatedDataFrame", fd), 
                        phenoData = new("NAnnotatedDataFrame", 
                                        data.frame(sampleNames = "1")), 
                        experimentData = new("MIAPE", 
                                             instrumentManufacturer = "a", 
                                             instrumentModel = "a", 
                                             ionSource = "a", 
                                             analyser = "a", 
                                             detectorType = "a"))
    pkAreaL[[i]] <- pkArea[is.na(pkGrpVal[, i]), , drop = FALSE]
  }
  rm(pkGrpVal)
  rm(pkArea)
  rm(min_fdata)
  message(" OK\nStart integrating peak areas from original files")
  ph <- processHistory(object, type = xcms:::.PROCSTEP.PEAK.DETECTION)
  findPeakMethod <- "unknown"
  mzCenterFun <- "wMean"
  if (length(ph)) {
    if (is(ph[[1]], "XProcessHistory")) {
      prm <- ph[[1]]@param
      findPeakMethod <- .param2string(prm)
      if (.hasSlot(prm, "mzCenterFun")) 
        mzCenterFun <- prm@mzCenterFun
    }
  }
  cp_colnames <- colnames(chromPeaks(object))
  mzCenterFun <- paste("mzCenter", gsub(mzCenterFun, 
                                        pattern = "mzCenter.", replacement = "", 
                                        fixed = TRUE), sep = ".")
  if (findPeakMethod == "MSW") {
    rts <- rtime(object, bySample = TRUE)
    if (any(lengths(rts) > 1)) 
      stop("The data is supposed to be direct injection data, ", 
           "but I got files with more than one spectrum/", 
           "retention time!")
    res <- bpmapply(FUN = xcms:::.getMSWPeakData, objectL, pkAreaL, 
                    as.list(1:length(objectL)), MoreArgs = list(cn = cp_colnames), 
                    BPPARAM = BPPARAM, SIMPLIFY = FALSE)
  } else if (findPeakMethod == "matchedFilter") {
    res <- bpmapply(FUN = xcms:::.getChromPeakData_matchedFilter, 
                    objectL, pkAreaL, as.list(1:length(objectL)), 
                    MoreArgs = list(cn = cp_colnames, param = prm, 
                                    msLevel = msLevel), 
                    BPPARAM = BPPARAM, 
                    SIMPLIFY = FALSE)
  } else {
    res <- bpmapply(FUN = xcms:::.getChromPeakData, objectL, 
                    pkAreaL, as.list(1:length(objectL)), 
                    MoreArgs = list(cn = cp_colnames, mzCenterFun = mzCenterFun, 
                                    msLevel = msLevel), 
                    BPPARAM = BPPARAM, SIMPLIFY = FALSE)
  }
  rm(objectL)
  res <- do.call(rbind, res)
  res <- cbind(res, group_idx = unlist(lapply(pkAreaL, function(z){
    z[, "group_idx"]}), use.names = FALSE))
  res <- res[!is.na(res[, "into"]), , drop = FALSE]
  if (nrow(res) == 0) {
    warning("Could not integrate any signal for the missing ", 
            "peaks! Consider increasing 'expandMz' and 'expandRt'.")
    return(object)
  }
  rm(pkAreaL)
  gc()
  newFd <- new("MsFeatureData")
  newFd@.xData <- xcms:::.copy_env(object@msFeatureData)
  object@msFeatureData <- new("MsFeatureData")
  incr <- nrow(chromPeaks(newFd))
  for (i in unique(res[, "group_idx"])) {
    fdef$peakidx[[i]] <- c(fdef$peakidx[[i]], (which(res[, "group_idx"] == i) + incr))
  }
  fdef <- rbind(fdef, featureDefinitions(newFd)[featureDefinitions(newFd)$ms_level != 
                                                  msLevel, , drop = FALSE])
  if (!any(colnames(fdef) == "ms_level")) {
    fdef$ms_level <- 1L
  } else {fdef <- fdef[order(fdef$ms_level), ]}
  reg_extract <- rownames(chromPeaks(newFd))
  maxId <- max(as.numeric(regmatches(reg_extract,regexpr("\\d+$", reg_extract))))
  if(is.na(maxId)){
    message("Grepping for max compound number dun fuxed, plz fix")
    browser()
  }
  if (maxId < 1)
    stop("chromPeaks matrix lacks rownames; please update ", 
         "'object' with the 'updateObject' function.")
  toId <- maxId + nrow(res)
  rownames(res) <- sprintf(paste0("CP", "%0", ceiling(log10(toId + 1L)), "d"), 
                           (maxId + 1L):toId)
  chromPeaks(newFd) <- rbind(chromPeaks(newFd), res[, -ncol(res)])
  cpd <- chromPeakData(newFd)[rep(1L, nrow(res)), , drop = FALSE]
  cpd[, ] <- NA
  cpd$ms_level <- as.integer(msLevel)
  cpd$is_filled <- TRUE
  if (!any(colnames(chromPeakData(newFd)) == "is_filled")) 
    chromPeakData(newFd)$is_filled <- FALSE
  chromPeakData(newFd) <- rbind(chromPeakData(newFd), cpd)
  rownames(chromPeakData(newFd)) <- rownames(chromPeaks(newFd))
  featureDefinitions(newFd) <- fdef
  lockEnvironment(newFd, bindings = TRUE)
  object@msFeatureData <- newFd
  ph <- xcms:::XProcessHistory(param = param, date. = startDate, 
                               type. = xcms:::.PROCSTEP.PEAK.FILLING, 
                               fileIndex = 1:length(fileNames(object)), 
                               msLevel = msLevel)
  object <- xcms:::addProcessHistory(object, ph)
  object
}
formula2elements <- function(formula_vec){
  split_formulas <- formula_vec %>%
    gregexpr(pattern = "[A-Z][a-z]*[0-9]*") %>%
    regmatches(x = formula_vec)
  elements_in <- lapply(split_formulas, gsub, pattern = "[0-9]", replacement = "")
  element_counts <- lapply(split_formulas, gsub, pattern = "[A-z]", replacement = "")
  element_counts <- lapply(element_counts, function(x){
    x[!nchar(x)]<-1
    return(as.numeric(x))
  })
  mapply(`names<-`, element_counts, elements_in, SIMPLIFY = FALSE)
}
grabSingleFileMS2 <- function(filename){
  msdata <- mzR::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  ms2rows <- seq_len(nrow(fullhd))[fullhd$msLevel>1]
  spectra_list <- lapply(ms2rows, function(x){
    rtime <- fullhd[x, "retentionTime"]
    premz <- fullhd[x, "precursorMZ"]
    fragments <- mzR::peaks(msdata, x)
    return(cbind(rtime, premz, fragments))
  })
  all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)), 
                        c("rt", "premz", "fragmz", "int"))
  return(all_data)
}
mgf_maker <- function(feature_msdata, ms1, ms2, output_file){
  ms1_to_write <- paste("BEGIN IONS",
                         paste0("PEPMASS=", feature_msdata$mzmed),
                         "MSLEVEL=1",
                         "CHARGE=1+",
                         apply(ms1, 1, paste, collapse=" "),
                         "END IONS",
                         "", sep = "\n")
  
  if(!nrow(ms2)){
    outtext <- ms1_to_write
  } else {
    ms2_to_write <- lapply(unique(ms2$voltage), function(volt_i){
      paste0(c("BEGIN IONS",
             paste0("PEPMASS=", feature_msdata$mzmed),
             "MSLEVEL=2",
             "CHARGE=1+",
             paste0(apply(ms2[ms2$voltage==volt_i, c("fragmz", "int")], 
                   1, paste, collapse=" "), collapse="\n"),
             "END IONS",
             ""), collapse="\n")
    })
    
    outtext <- paste0(ms1_to_write, "\n", paste0(unlist(ms2_to_write), collapse = "\n"))
  }
  writeLines(outtext, con = output_file)
}
isocheck <- function(feature_num, final_peaks=final_peaks, printplot=FALSE){
  ft_isodata <- final_peaks %>% filter(feature==feature_num) %>%
    select(M_area, filename, C13_area, N15_area, O18_area, X2C13_area, S34_area, S33_area) %>%
    pivot_longer(cols = starts_with(c("C", "X", "N", "O", "S")), names_to="isotope") %>%
    filter(value!=0) %>%
    mutate(isotope=gsub(pattern = "_area", replacement = "", x = .$isotope))
  
  if(all(is.na(ft_isodata))){
    return(data.frame(
      isotope=c("C13", "N15", "O18", "X2C13", "S33", "S34"),
      count=NA
      ) %>% `colnames<-`(c("isotope", feature_num)))
  }
  
  count_ests <- split(ft_isodata, ft_isodata$isotope) %>%
    map(lm, formula=value~M_area) %>%
    map(summary) %>%
    imap(function(x, y){
      iso_abundance <- iso_abundance_table %>%
        filter(isotope==y)
      
      rsquared <- x$r.squared
      if(is.na(rsquared))return(NA)
      if(rsquared<0.99|rsquared==1)return(NA)
      coefs <- x$coefficients
      norm_factor <- sapply(0:10, dbinom, x=0, 
                            prob=iso_abundance$abundance)
      pred_values <- sapply(0:10, dbinom, x=iso_abundance$n_atoms, 
                            prob=iso_abundance$abundance)
      est_slope <- coefs["M_area", "Estimate"]
      (0:10)[which.min(abs(pred_values/norm_factor-est_slope))]
    })
  count_ests <- data.frame(names(count_ests), unlist(count_ests)) %>%
    `rownames<-`(NULL) %>% `colnames<-`(c("isotope", feature_num)) %>%
    filter(ft_isodata %>% 
             group_by(isotope) %>% 
             summarize(count=n()) %>% 
             pull(count) %>% 
             `>=`(5))
  
  if(printplot){
    gp <- ggplot(ft_isodata, aes(x=M_area, y=value)) + 
      geom_point() + 
      geom_smooth(method = "lm") +
      facet_wrap(~isotope, scales = "free") +
      ggtitle(feature_num)
    print(gp)
  }
  
  return(count_ests)
}
rdisop_check <- function(feature_num, feature_data, database_formulae){
  feature_msdata <- feature_data[feature_data$feature==feature_num, ]
  ms1 <- rbind(c(feature_msdata$mzmed, feature_msdata$avgarea),
               c(feature_msdata$mzmed+1.003355, feature_msdata$C13),
               c(feature_msdata$mzmed+1.003355*2, feature_msdata$X2C13),
               c(feature_msdata$mzmed+1.995796, feature_msdata$S34),
               c(feature_msdata$mzmed+0.997035, feature_msdata$N15),
               c(feature_msdata$mzmed+2.004244, feature_msdata$O18))
  ms1 <- ms1[ms1[,2]!=0, , drop=FALSE]
  # ms1[,1] <- ms1[,1]-1.007276
  ms1[,1] <- ms1[,1] # Try this and see what changes? Assumes Rdisop has wrong e- mass
  rdoutput <- Rdisop::decomposeIsotopes(masses = ms1[,1], intensities = ms1[,2], 
                                        ppm = ifelse(ms1[1,1]<200, 5*200/ms1[1,1], 5), 
                                        maxisotopes = 2)
  if(is.null(rdoutput)){return(NA)}
  rd_df <- rdoutput %>%
    `[[<-`("isotopes", NULL) %>%
    do.call(what = cbind) %>% 
    as.data.frame(stringsAsFactors=FALSE) %>%
    #filter(valid=="Valid"&DBE>=-1) %>%
    filter(formula%chin%database_formulae)
  if(!nrow(rd_df)){return(NA)}
  rd_df %>% #slice(1) %>%
    pull(formula) %>%
    #unlist()
    list()
}
