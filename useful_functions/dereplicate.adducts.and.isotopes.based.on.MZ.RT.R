## The function spits out a data frame listing pairs of mass features that are adducts/isotopologues
## the input needs the columns c(MassFeature, Average.Mz, Average.Rt.min.) 

## The function spits out a data frame listing pairs of mass features that are adducts/isotopologues
## the input dataframe needs the columns c(MassFeature, Average.Mz, Average.Rt.min.) 
## you also need to specify ionization mode, MZ threshold and RT threshold to be 
## considered candidates for adducts/isotopologues

## de-replicate adducts by pariwise RT MZ comparisons----------------
dereplicate.adducts = function(data, MZthreshold = 0.0025, RTthreshold = 0.2, ionization = "pos") {
  require(ggplot2)
  require(tidyr)
  require(stringr)
  require(readr)
  require(dplyr)
  require(vegan)
  require(tidyverse)
  require(gtools)
  require(cluster)
  require(factoextra)
  require(viridis)
  
  data.1 <- data %>%
    dplyr::select(MassFeature, Average.Mz, Average.Rt.min.) %>%
    unique()
  # add the key column
  controlTable <- data.frame(combinations(n=length(data.1$MassFeature),
                                          r=2,
                                          v=c(data.1$MassFeature),
                                          repeats.allowed = FALSE))
  colnames(controlTable) <- c("key1", "key2")
  controlTable <- cbind(
    data.frame(pair_key = paste(controlTable[[1]], controlTable[[2]]),
               stringsAsFactors = FALSE),
    controlTable) 
  if(ionization=="pos"){
    possible.mass.offsets <- c(
      17.02655, # difference between M+H and M+NH4
      21.98194, # difference between M+H and M+Na
      37.95588, # difference between M+H and M+K
      1.003354, #13C
      2*1.003354, # 2 x 13C
      1.006277, #2H
      0.997036, #15N
      2.004245, #18O
      18.010565, # plus H2O
      18.010565*2, # plus 2xH2O
      1.997306, #32S vs 34S
      0.998626, #32S vs33S
      41.02655, ## plus ACN
      32.02621, #plus MeOH
      63.00844 ## M+ACN+Na
    )
  } else {
    possible.mass.offsets <- c(
      19.96739, # difference between M-H and M+Na-2H
      35.97668, # difference between M-H and M+Cl
      37.95643, # difference between M-H and M+K-2H
      1.003354, #13C
      2*1.003354, # 2 x 13C
      1.006277, #2H
      0.997036, #15N
      2.004245, #18O
      18.01111, # plus/minus H2O
      18.010565*2, # plus/minus 2xH2O
      1.997306, #32S vs 34S
      0.998626, #32S vs33S
      79.92616 # difference between M-H and M+Br
    )
  }
  
  bad.apple.contestants <- full_join(controlTable, data.1, by = c("key1"="MassFeature")) %>%
    dplyr::rename(Average.Mz.1 = Average.Mz,
                  Average.Rt.min.1 = Average.Rt.min.) %>%
    full_join(., data.1, by = c("key2"="MassFeature")) %>%
    dplyr::rename(Average.Mz.2 = Average.Mz,
                  Average.Rt.min.2 = Average.Rt.min.) %>%
    filter(!is.na(pair_key)) %>%
    mutate(Mass.Diff = abs(Average.Mz.1-Average.Mz.2),
           RT.Diff = abs(Average.Rt.min.1-Average.Rt.min.2)) %>%
    filter(RT.Diff < RTthreshold) 
  MZ.test <- bad.apple.contestants$Mass.Diff
  names(MZ.test) <- bad.apple.contestants$pair_key
  MZ.test.results <- data.frame(t(sapply(MZ.test, function(x) x- possible.mass.offsets)))
  MZ.test.results <- MZ.test.results %>%
    mutate(pair_key = rownames(MZ.test.results))
  bad.apple.contestants.MZ.filter <- bad.apple.contestants %>%
    full_join(MZ.test.results) %>%
    filter(!is.na(key1)) %>%
    gather(MZ.diffs.to.compare, Diff.Val, -pair_key, -key1,-key2,-Average.Mz.1,
           -Average.Rt.min.1, -Average.Mz.2, -Average.Rt.min.2, -Mass.Diff, -RT.Diff) %>%
    mutate(Diff.Val = abs(Diff.Val)) %>%
    filter(Diff.Val < MZthreshold)
  bad.apple.contestants.refined <- bad.apple.contestants %>%
    filter(pair_key %in% bad.apple.contestants.MZ.filter$pair_key) %>%
    mutate()
  bad.apple.contestants.refined
}

testing <- dereplicate.adducts(data = "data/Area_HILICNeg_EddyTransect.csv", MZthreshold = 0.0025, RTthreshold = 0.2, ionization = "pos")