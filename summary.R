library(knitr)
# Function to extract desired parameters from a saved file
extract_params <- function(file_path) {
  load(file_path)
  
  # Extract desired parameters
  params <- list(
    f1_jive = results.save$results['jive', "f1"],
    f1_slide = results.save$results['slide', "f1"],
    f1_ajive = results.save$results['ajive', "f1"],
    f1_dcca = results.save$results['dcca', "f1"],
    f1_unifac = results.save$results['unifac', "f1"],
    f1_proposed = results.save$results['proposed', "f1"],
    f1_proposed_subsampling = results.save$results['proposed_subsampling', "f1"],
    phi.max = results.save$phi_max,
    SNR1 = results.save$SNR1,
    SNR2 = results.save$SNR2,
    no_joint = results.save$no_joint,
    no_indiv = results.save$no_indiv
  )
  
  return(params)
}

# Get paths of all saved files starting with "demo2_"
file_paths <- list.files(path = './results', pattern = "^demo2_.*\\.RData$", full.names = TRUE)

# Extract parameters from each file
all_params <- lapply(file_paths, extract_params)

# Convert list of lists to data frame
df <- do.call(rbind, lapply(all_params, data.frame, row.names = NULL))

# Organize the data frame
df <- df[, c("f1_jive",
             "f1_slide",
             "f1_ajive", 
             "f1_dcca",
             "f1_unifac",
             "f1_proposed", 
             "f1_proposed_subsampling",
             "phi.max", 
             "SNR1", 
             "SNR2", 
             "no_joint", 
             "no_indiv")]

# Print the data frame
df_short <- df[1:9, c("f1_jive", "f1_slide", "f1_ajive", "f1_dcca", "f1_unifac", "f1_proposed", "f1_proposed_subsampling")]
df_short <- cbind(expand.grid(c('high SNR', 'medium SNR', 'low SNR'), c('Orthogonal', 'Slightly aligned', 'Aligned')), df_short)
# round to 3 decimal places
df_short[, c('f1_jive', 'f1_slide', 'f1_ajive', 'f1_dcca', 'f1_unifac', 'f1_proposed', 'f1_proposed_subsampling')] <- round(df_short[, c('f1_jive', 'f1_slide', 'f1_ajive', 'f1_dcca', 'f1_unifac', 'f1_proposed', 'f1_proposed_subsampling')], 3)
kable(df_short[, c('Var2', 'Var1', 'f1_jive', 'f1_slide', 'f1_ajive', 'f1_dcca', 'f1_unifac', 'f1_proposed', 'f1_proposed_subsampling')], 'latex')
