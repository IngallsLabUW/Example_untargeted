
# Custom function that tries to resolve discrepancies between various metrics
# Basically a bunch of if/else statements
# Highest priority: isotope validated (isotope_choice)
# Next highest: if two peaks are found and two standards should be found, use RT
#               ordering to match them up (match_choice)
# Next highest: if peaks are dramatically larger in the correct mix and smaller
#               in the other mix, assume it's the larger one (mix_choice)
#               sometimes produces multiple, thus later area and RT matching
# Next highest: if one peak is larger than the others, assume it's the
#               standard (area_choice)
# Lowest: RT matching based on the number in the standards list
#
# Returns a character vector, usually a single peak but can be multiple or none
# concatenated with a semicolon
stan_guesser <- function(isotope_choice, mix_choice, match_choice, area_choice, rt_choice){
  if(all(is.na(c(isotope_choice, mix_choice, match_choice, area_choice, rt_choice)))){
    return("No peaks found")
  }
  if(!is.na(isotope_choice)){
    return(isotope_choice)
  }
  if(!is.na(match_choice)) {
    return(match_choice)
  }
  if(is.na(mix_choice)){
    mix_options <- NA
  } else if(!nchar(mix_choice)){
    mix_options <- NA
  } else {
    mix_options <- unlist(strsplit(mix_choice, split = "; "))
  }
  if(length(mix_options)==1&!all(is.na(mix_options))){
    return(mix_options)
  }
  
  if(length(mix_options)>1&!all(is.na(mix_options))){
    if(rt_choice==area_choice){
      return(rt_choice)
    }
    if(rt_choice%in%mix_options&!area_choice%in%mix_options){
      return(rt_choice)
    }
    if(area_choice%in%mix_options&!rt_choice%in%mix_options){
      return(area_choice)
    }
    if(area_choice%in%mix_options&rt_choice%in%mix_options){
      return(paste(c(area_choice, rt_choice), collapse = "; "))
    }
  }
  return(rt_choice)
}

# For each standard
# Grab all the features that it could be based on 5ppm m/z window
# If there aren't any features, return NA for everything

# If there are features, calculate 5 metrics for each one:
# Isotope choice:
#   If the compound has an isotope, use the peak closest to the isotope peak in RT
#     (i.e. another compound exists with the same name + C13, N15, or O18)
#   If the compound IS an isotope (i.e. has C13, N15, or O18 in the name)
#     find the two peaks closest in RT and assume those are the isotopologues
#   Otherwise, NA
# Match choice:
#   If there are two standards expected in a mix and two features found,
#   assume that the two found peaks are the standards and use rank-order RT matching
#   If more or fewer peaks found, NA
# Mix choice:
#   Generally, which peaks are bigger in the mix they've been added to?
#   Calculate average z-statistic between the mixes after controlling for H2O vs Matrix
#   Return all feature numbers with a z-score above 10 (arbitrarily chosen)
# Area choice:
#   Which peak has the largest area?
# RT choice:
#   Which peak is closest in RT to the value in the standards list?

# Then use function above to resolve discrepancies based on the priority hierarchy
# Manual assignments for known annoyances at the bottom

stan_annotations <- given_stans %>%
  split(.$compound_name) %>%
  pblapply(function(stan_data){
    # dput(stan_data)
    possible_stan_features <- feature_data %>%
      filter(mzmed%between%pmppm(as.numeric(stan_data["mz"]), 6)) %>%
      arrange(rtmed)
    if(!nrow(possible_stan_features)){
      return(rep(NA, 5))
    }
    
    # If the compound has an isotopologue, use RT matching
    # If the compound IS an isotopologue, same deal
    isotopologue <- stan_data$compound_name %>%
      paste0(", [15N\\d|13C\\d|2H\\d]") %>%
      grep(x = given_stans$compound_name, value = TRUE)
    if(length(isotopologue)==1){
      isotopo_data <- given_stans %>% 
        filter(compound_name==isotopologue)
      isotopo_features <- feature_data %>%
        filter(mzmed%between%pmppm(as.numeric(isotopo_data$mz), 6))
      if(nrow(isotopo_features)==1){
        isotope_choice <- possible_stan_features %>% 
          mutate(rtdiff=abs(rtmed-isotopo_features$rtmed)) %>%
          arrange(rtdiff) %>%
          slice(1) %>%
          pull(feature)
      } else {
        isotope_choice <- NA
      }
    } else if(grepl(pattern = ", 15N|13C|2H", x = stan_data$compound_name)) {
      isotopologue <- gsub(", 15N.*|, 13C.*|, 2H.*", "", stan_data$compound_name)
      isotopo_data <- given_stans %>% 
        filter(compound_name==isotopologue|
                 compound_name==gsub("^D", "", isotopologue)|
                 compound_name==gsub("^L", "", isotopologue))
      isotopo_features <- feature_data %>%
        filter(mzmed%between%pmppm(as.numeric(isotopo_data$mz), 5))
      isotope_choice <- possible_stan_features %>% 
        left_join(isotopo_features, by=character()) %>%
        select(-starts_with("mzmed")) %>%
        mutate(rtdiff=abs(rtmed.x-rtmed.y)) %>%
        arrange(rtdiff) %>%
        slice(1) %>%
        pull(feature.x)
    } else {
      isotope_choice <- NA
    }
    
    # If there's a peak that differs between the mixes
    if(is.na(stan_data$mix)){
      mix_choice <- NA
    } else {
      stan_peaks <- filled_peaks %>%
        filter(feature%in%possible_stan_features$feature) %>%
        filter(str_detect(filename, "Mix")) %>%
        filter(grepl("Std", filename)) %>%
        filter(!grepl("H2OinMatrix", filename)) %>%
        mutate(date_run=strptime(str_extract(filename, "^\\d+"), format = "%y%m%d")) %>%
        filter(date_run>=stan_data$date_added) %>%
        select(feature, mz, rt, M_area, filename) %>%
        mutate(correct_mix=grepl(filename, pattern = stan_data["mix"])) %>%
        mutate(stan_type=str_extract(filename, "InH2O|InMatrix"))
      
      if(!nrow(stan_peaks)){
        mix_choice <- NA
      } else {
        mix_peaks <- stan_peaks %>% 
          group_by(feature, correct_mix, stan_type) %>%
          summarise(avgarea=mean(M_area), sdarea=sd(M_area)) %>%
          right_join(expand.grid(
            unique(.$feature),
            unique(.$correct_mix),
            unique(.$stan_type)
          ) %>% `names<-`(c("feature", "correct_mix", "stan_type")),
          by=c("feature", "correct_mix", "stan_type")) %>%
          mutate(avgarea=ifelse(is.na(avgarea), 0, avgarea)) %>%
          mutate(sdarea=ifelse(is.na(sdarea), 0, sdarea))
        
        if(length(unique(mix_peaks$feature))==1){
          mix_choice <- unique(mix_peaks$feature)
        } else {
          mix_choice <- mix_peaks %>%
            split(interaction(.$feature, .$stan_type)) %>%
            lapply(function(v){
              if(nrow(v)==1){
                if(v$correct_mix==TRUE){
                  data.frame(feature=unique(v$feature), 
                             stan_type=unique(v$stan_type),
                             diff_degree=Inf)
                } else {
                  data.frame(feature=unique(v$feature), 
                             stan_type=unique(v$stan_type),
                             diff_degree=-Inf)
                }
              } else {
                diff <- (v$avgarea[v$correct_mix] - v$avgarea[!v$correct_mix])/mean(v$sdarea)
                data.frame(feature=unique(v$feature), 
                           stan_type=unique(v$stan_type),
                           diff_degree=diff)
              }
            }) %>%
            do.call(what = rbind) %>% `rownames<-`(NULL) %>%
            group_by(feature) %>%
            summarise(correct_mix_peak=mean(diff_degree)) %>%
            filter(correct_mix_peak>1) %>%
            pull(feature) %>%
            paste(collapse = "; ")
        }
      }
    }
    if(!is.na(mix_choice))if(!nchar(mix_choice))mix_choice<-NA

    # If there's one peak much closer in RT to expected than the others
    expected_rt <- stan_data$rt*60
    rt_choice <- possible_stan_features %>% 
      mutate(rtdiff=abs(rtmed-expected_rt)) %>%
      arrange(rtdiff) %>%
      slice(1) %>%
      pull(feature)
    
    # If there's same number of features as expected peaks, assume 1:1 and order by RT
    possible_other_stans <- given_stans %>%
      filter(mz%between%pmppm(stan_data$mz, 5)) %>%
      arrange(rt)
    if(nrow(possible_stan_features)==nrow(possible_other_stans)){
      match_choice <- possible_stan_features %>%
        arrange(rtmed) %>%
        cbind(possible_other_stans) %>%
        select(feature, compound_name) %>%
        filter(compound_name==stan_data$compound_name) %>%
        pull(feature)
    } else {
      match_choice <- NA
    }
    
    stan_peaks <- filled_peaks %>%
      filter(feature%in%possible_stan_features$feature) %>%
      select(feature, mz, rt, M_area, filename) %>%
      filter(grepl("Std", filename))
    
    if(nrow(stan_peaks)==0){
      area_choice <- NA
    } else if(nrow(possible_stan_features)==nrow(possible_other_stans)&
              nrow(possible_stan_features)>1) {
      area_choice <- NA
    } else {
      area_choice <- stan_peaks %>% 
        group_by(feature) %>%
        summarise(avgarea=mean(M_area)) %>%
        arrange(desc(avgarea)) %>%
        slice(1) %>%
        pull(feature)
    }
    
    feature <- stan_guesser(isotope_choice, mix_choice, match_choice, area_choice, rt_choice)
    
    return(data.frame(
      feature=feature,
      isotope_validated=isotope_choice,
      rt_matchup=match_choice,
      mix_matched=mix_choice,
      closer_rt=rt_choice,
      area_choice=area_choice
    ))
  }) %>%
  do.call(what=rbind) %>%
  mutate(compound_name=rownames(.)) %>%
  select(compound_name, everything())


# Check on internal standard assignments
stan_annotations %>%
  filter(!is.na(.$isotope_validated)&isotope_validated!="Manual")


# Check on features that are dual-assigned - these should be resolved by manual assignment
stan_annotations[(duplicated(stan_annotations$feature, fromLast=TRUE)|
                   duplicated(stan_annotations$feature))&
                   !is.na(stan_annotations$feature), ] %>%
  split(.$feature)


# Manual assignments
if(cruise=="falkor"){
  glutamate_data <- given_stans %>% filter(compound_name=="L-Glutamic acid")
  glutamate <- feature_data %>%
    filter(mzmed%between%pmppm(glutamate_data$mz, 5)) %>%
    arrange(desc(avgarea)) %>%
    slice(1)
  stan_annotations[stan_annotations$compound_name=="L-Glutamic acid",] <- 
    c("L-Glutamic acid", glutamate$feature, rep("Manual", ncol(stan_annotations)-2))
  stan_annotations[stan_annotations$compound_name=="beta-Glutamic acid",] <- 
    c("beta-Glutamic acid", rep(NA, ncol(stan_annotations)-1))
  stan_annotations[stan_annotations$compound_name=="O-Acetyl-L-serine",] <- 
    c("O-Acetyl-L-serine", rep(NA, ncol(stan_annotations)-1))
  
  threo_data <- given_stans %>% filter(compound_name=="L-Threonine")
  threo <- feature_data %>%
    filter(mzmed%between%pmppm(threo_data$mz, 5)) %>%
    arrange(desc(avgarea)) %>%
    slice(1)
  stan_annotations[stan_annotations$compound_name=="L-Threonine",] <- 
    c("L-Threonine", threo$feature, rep("Manual", ncol(stan_annotations)-2))
  stan_annotations[stan_annotations$compound_name=="L-Homoserine",] <- 
    c("L-Homoserine", rep(NA, ncol(stan_annotations)-1))
  
  ach_data <- given_stans %>% filter(compound_name=="Acetylcholine")
  stan_annotations[stan_annotations$compound_name=="Acetylcholine",] <-
    c("Acetylcholine", "FT0617", rep("Manual", ncol(stan_annotations)-2))
} else if(cruise=="mesotransect") {
  # Threonine is the left peak, homoserine is the right one
  # According to standards list RT column
  # However, the left peak seems to have not been detected??
  threo_data <- given_stans %>% filter(compound_name=="L-Threonine")
  threo <- feature_data %>%
    filter(mzmed%between%pmppm(threo_data$mz, 5)) %>%
    arrange(desc(avgarea)) %>%
    slice(1)
  stan_annotations[stan_annotations$compound_name=="L-Threonine",] <- 
    c("L-Threonine", threo$feature, rep("Manual", ncol(stan_annotations)-2))
  stan_annotations[stan_annotations$compound_name=="L-Homoserine",] <- 
    c("L-Homoserine", rep(NA, ncol(stan_annotations)-1))
  
  # We don't see isoleucine at all, just leucine (based on 15N isoleucine)
  stan_annotations[stan_annotations$compound_name=="L-Isoleucine",] <- 
    c("L-Isoleucine", rep(NA, ncol(stan_annotations)-1))
  
  # Same with valine - completely swallowed by betaine
  stan_annotations[stan_annotations$compound_name=="L-Valine",] <- 
    c("L-Valine", rep(NA, ncol(stan_annotations)-1))
  
  stan_annotations[stan_annotations$compound_name=="O-Acetyl-L-serine",] <- 
    c("O-Acetyl-L-serine", rep(NA, ncol(stan_annotations)-1))
  
} else if(cruise=="mesocenter"){
  threo_data <- given_stans %>% filter(compound_name=="L-Threonine")
  threo <- feature_data %>%
    filter(mzmed%between%pmppm(threo_data$mz, 5)) %>%
    arrange(desc(avgarea)) %>%
    slice(1)
  stan_annotations[stan_annotations$compound_name=="L-Threonine",] <- 
    c("L-Threonine", threo$feature, rep("Manual", ncol(stan_annotations)-2))
  stan_annotations[stan_annotations$compound_name=="L-Homoserine",] <- 
    c("L-Homoserine", rep(NA, ncol(stan_annotations)-1))
  
  # We don't see isoleucine at all, just leucine (based on 15N isoleucine)
  stan_annotations[stan_annotations$compound_name=="L-Isoleucine",] <- 
    c("L-Isoleucine", rep(NA, ncol(stan_annotations)-1))

  # Same with valine - completely swallowed by betaine
  stan_annotations[stan_annotations$compound_name=="L-Valine",] <- 
    c("L-Valine", rep(NA, ncol(stan_annotations)-1))

  stan_annotations[stan_annotations$compound_name=="O-Acetyl-L-serine",] <- 
    c("O-Acetyl-L-serine", rep(NA, ncol(stan_annotations)-1))
}

stan_annotations[(duplicated(stan_annotations$feature, fromLast=TRUE)|
                    duplicated(stan_annotations$feature))&
                   !is.na(stan_annotations$feature), ] %>%
  split(.$feature) %>%
  print()
