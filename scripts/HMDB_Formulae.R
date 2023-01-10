library(xml2)

# Gotten by right clicking "Download" button and saying "Copy link address"
hmdb_metab_url <- "https://hmdb.ca/system/downloads/4.0/hmdb_metabolites.zip"
output_file <- "~/../Downloads/hmdb_metabolites.zip"

# File is ~600MB while compressed
download.file(hmdb_metab_url, destfile = output_file)

# Decompress to ~4GB
utils::unzip(zipfile = output_file, exdir = "~/../Downloads")

# Read in XML document (huge, requires hella memory)
v <- read_xml("~/../Downloads/hmdb_metabolites.xml")

# Determine node structure
xml_children(v)
# Looks like each metabolite has a single entry
# Look at structure within the first node
xml_children(v)[[1]]
# Looks like most interesting bits of data have their own node and
# can be extracted by name

# Get all the formulae
# d1: is a default namespace
# We want all the text within the <chemical_formula> tags
formulae <- xml_text(xml_find_all(v, "//d1:chemical_formula"))
# Get all the names
chem_names <- xml_text(xml_find_all(v, "//d1:name"))
# Can repeat for other nodes

# Convert to output df
output_df <- data.frame(chem_names, formulae)

# Save as compressed object to send to myself via Slack lol
saveRDS(formulae, file = "~/../Downloads/hmdb_formulae.rds")
Footer
© 2023 GitHub, Inc.
Footer navigation
Terms
Privacy
Security
Status
