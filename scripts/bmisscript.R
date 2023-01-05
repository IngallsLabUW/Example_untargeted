# Script to find the B-MIS for each standard in the untargeted data set
# Called by Control.Rmd

# all_peaks, given_stans, polarity, pmppm all defined in Control.Rmd


# Grab the internal standards and clean up a little
found_stans <- given_stans %>%
  filter(compound_type=="Internal Standard") %>%
  select(compound_name) %>%
  left_join(stan_annotations)


# Manually correct duplicate peaks
if(length(found_stans$compound_name)!=length(unique(found_stans$compound_name))){
  stop("Fix yo internal standard annotations!")
}


# Perform BMIS
pooled <- filled_peaks %>%
  filter(str_detect(filename, "_Poo_"))

pooled_IS <- pooled %>%
  filter(feature%in%found_stans$feature) %>%
  left_join(found_stans[,c("feature", "compound_name")]) %>%
  select(feature, filename, M_area, compound_name) %>%
  rename_with(~paste0("IS_", .), -filename) %>%
  add_row(filename=unique(.$filename), IS_compound_name="None", IS_M_area=1)

initial_CVs <- pooled %>%
  mutate(M_area=ifelse(str_detect(filename, "Half"), M_area*2, M_area)) %>%
  group_by(feature) %>%
  summarize(init_CV=cv(M_area))

standardized_CVs <- full_join(pooled, pooled_IS, by=c("filename")) %>%
  group_by(feature, IS_compound_name) %>%
  summarize(norm_CV=cv(M_area/IS_M_area), .groups="drop")

chosen_BMIS <- standardized_CVs %>%
  left_join(initial_CVs) %>%
  group_by(feature) %>%
  mutate(MIS_improvement=(init_CV - norm_CV)/init_CV) %>% 
  filter(MIS_improvement>min_improvement|IS_compound_name=="None") %>%
  arrange(desc(MIS_improvement)) %>%
  group_by(feature) %>%
  slice(1) %>%
  select(feature, BMIS=IS_compound_name)

IS_peaks <- found_stans %>%
  select(compound_name, feature) %>%
  left_join(filled_peaks, by="feature") %>%
  add_row(feature="FT0000", compound_name="None", filename=unique(.$filename), 
          M_area=ifelse(str_detect(filename, "Half"), 0.5, 1))
