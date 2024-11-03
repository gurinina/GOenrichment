## code to prepare datasets for GOenrichment package


# Load necessary libraries
library(usethis)
library(fgsea)

# Load the raw data
hGOBP.gmt <- gmtPathways("data-raw/hGOBP.gmt")

# Save processed data in the `data` directory
usethis::use_data(hGOBP.gmt, overwrite = TRUE)
