# Function to extract desired parameters from a saved file
extract_params <- function(file_path) {
  load(file_path)
  
  # Extract desired parameters
  params <- list(
    f1_ajive = results.save$results['ajive', "f1"],
    f1_proposed = results.save$results['proposed', "f1"],
    f1_proposed_subsampling1 = results.save$results['proposed_subsampling1', "f1"],
    f1_proposed_subsampling2 = results.save$results['proposed_subsampling2', "f1"],
    phi.max = results.save$phi_max,
    SNR1 = results.save$SNR1,
    SNR2 = results.save$SNR2,
    no_joint = results.save$no_joint,
    no_indiv = results.save$no_indiv
  )
  
  return(params)
}

# Get paths of all saved files starting with "demo2_"
file_paths <- list.files(path = './results', pattern = "^demo3_.*\\.RData$", full.names = TRUE)

# Extract parameters from each file
all_params <- lapply(file_paths, extract_params)

# Convert list of lists to data frame
df <- do.call(rbind, lapply(all_params, data.frame, row.names = NULL))

# Organize the data frame
df <- df[, c("f1_ajive", "f1_proposed", "f1_proposed_subsampling1", "f1_proposed_subsampling2", "phi.max", "SNR1", "SNR2", "no_joint", "no_indiv")]

# Print the data frame
print(df)