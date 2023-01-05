
message("Running SIRIUS formulas...")
mzs <- feature_data$mzmed
sirius_cmd <- paste0('sirius --noCite decomp',
                     ' --mass=', paste(mzs, collapse = ","),
                     ' --ppm=5',
                     ' --elements=CHNOPS')
if(nchar(sirius_cmd)>8096){
  # Windows doesn't support commands longer than 8,096 characters
  n_chunks <- ceiling(nchar(sirius_cmd)/8096)
  mz_chunks <- split(mzs, cut(seq_along(mzs), n_chunks, labels = FALSE))
  sirius_formulas_i <- lapply(mz_chunks, function(mzs_i){
    sirius_cmd_i <- paste0('sirius --noCite decomp',
                         ' --mass=', paste(mzs_i, collapse = ","),
                         ' --ppm=5',
                         ' --elements=CHNOPS')
    sirius_output <- system(sirius_cmd_i, intern = TRUE)
    sirius_output <- sirius_output[grep(pattern = "m/z", sirius_output):length(sirius_output)]
    
    output_lengths <- sapply(sirius_output, nchar, USE.NAMES = FALSE)
    long_groups <- cumsum(output_lengths!=8095)-as.numeric(output_lengths!=8095)
    full_lines <- split(sirius_output, long_groups) %>%
      sapply(paste0, USE.NAMES = FALSE, collapse="")
    sirius_output <- grep(pattern = "^\\d", full_lines, value = TRUE)
    sirius_output <- gsub(sirius_output, pattern = "\\t$", replacement="\\\tNA")
    
    sirius_formulas <- read.table(text=sirius_output)
  })
  sirius_formulas <- do.call(sirius_formulas_i, what = "rbind")
} else {
  sirius_output <- system(sirius_cmd, intern = TRUE)
  sirius_output <- sirius_output[grep(pattern = "m/z", sirius_output):length(sirius_output)]
  
  output_lengths <- sapply(sirius_output, nchar, USE.NAMES = FALSE)
  long_groups <- cumsum(output_lengths!=8095)-as.numeric(output_lengths!=8095)
  full_lines <- split(sirius_output, long_groups) %>%
    sapply(paste0, USE.NAMES = FALSE, collapse="")
  sirius_output <- grep(pattern = "^\\d", full_lines, value = TRUE)
  sirius_output <- gsub(sirius_output, pattern = "\\t$", replacement="\\\tNA")
  
  sirius_formulas <- read.table(text=sirius_output)
}

sirius_formulas <- sirius_formulas[,2] %>%
  strsplit(split = ",")
names(sirius_formulas) <- feature_data$feature




# Run Rdisop ----
message("Running Rdisop...")
rdisop_formulas <- feature_data$feature %>%
  pbsapply(rdisop_check, feature_data = feature_data, 
           database_formulae = database_formulae)


# Check isotope matches ----
message("Running isotope checks...")
iso_abundance_table <- data.frame(
  isotope=c("C13", "N15", "O18", "S33", "S34", "X2C13"),
  abundance=c(0.011, 0.00368, 0.00205, 0.0075, 0.0421, 0.011),
  n_atoms=c(1,1,1,1,1,2)
)
iso_formulas <- feature_data$feature %>%
  pblapply(isocheck, final_peaks=peak_envelopes) %>%
  c(list(data.frame(isotope=c("C13", "N15", "O18", "X2C13", "S33", "S34"))), .) %>%
  purrr::reduce(.f=left_join, by="isotope") %>%
  t() %>% as.data.frame() %>% mutate(feature=rownames(.)) %>%
  `colnames<-`(slice(., 1) %>% `[`(-length(.)) %>% c("feature")) %>%
  slice(-1) %>% select(feature, everything())



# Assemble formulas, remove disagreeing ones and duplicates ----
inter_formulas <- sapply(names(rdisop_formulas), function(feature_num){
  intersect(rdisop_formulas[[feature_num]], sirius_formulas[[feature_num]])
}, simplify=FALSE)

isochecked_formulas <- lapply(names(inter_formulas), function(feature_num){
  isodata <- filter(iso_formulas, feature==feature_num)
  if(!length(inter_formulas[[feature_num]])){
    return(character(0))
  }
  formula_agreements <- inter_formulas[[feature_num]] %>%
    formula2elements() %>%
    sapply(function(elem_table){
      empty_elements <- c("C", "N", "O", "S")[!c("C", "N", "O", "S")%in%names(elem_table)]
      empty_elements <- `names<-`(numeric(length(empty_elements)), empty_elements)
      elem_table <- c(elem_table, empty_elements)
      c(
        C=elem_table[["C"]]-as.numeric(isodata[["C13"]]),
        N=elem_table[["N"]]-as.numeric(isodata[["N15"]]),
        O=elem_table[["O"]]-as.numeric(isodata[["O18"]]),
        S=elem_table[["S"]]-as.numeric(isodata[["S34"]])
      )
  }, simplify = FALSE) %>%
    sapply(sum, na.rm=TRUE) %>%
    sapply(abs)
  best_formula <- inter_formulas[[feature_num]][
    which(formula_agreements==min(formula_agreements))
  ]
  return(best_formula)
}) %>% `names<-`(names(inter_formulas))

feature_formulas <- data.frame(feature=names(isochecked_formulas),
           formula=sapply(isochecked_formulas, paste, collapse="; ")) %>%
  left_join(feature_data, by="feature") %>%
  select(feature, formula)
