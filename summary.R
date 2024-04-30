library(knitr)
# Function to extract desired parameters from a saved file

extract_params <- function(file_path) {
  load(file_path)
  
  out <- list()
  models <- names(results.save$results)
  for (model in models) {
    table <- results.save$results[[model]]
    table[['precision.joint']] <- 1 - table$fdr.Pjoint
    table[['precision.indiv1']] <- 1 - table$fdr.Pindiv1
    table[['precision.indiv2']] <- 1 - table$fdr.Pindiv2
    table[['avg.precision']] <- (table[['precision.joint']] + 
                                   table[['precision.indiv1']] + 
                                   table[['precision.indiv2']]) / 3
    table[['avg.tpr']] <- (table[['tpr.Pjoint']] + 
                             table[['tpr.Pindiv1']] + 
                             table[['tpr.Pindiv2']]) / 3
    table[['f1.avg']] <- 2 * table$avg.precision * table$avg.tpr / (table$avg.precision + table$avg.tpr)
    out[[paste0(model, '_avgF1')]] <- mean(table$f1.avg)
    out[[paste0(model, '_stdF1')]] <- sd(table$f1.avg)
  }
  
  # Extract desired parameters
  out[['phi_max']] <- results.save$phi_max
  out[['SNR1']] <- results.save$SNR1
  out[['SNR2']] <- results.save$SNR2
  out[['no_joint']] <- results.save$no_joint
  out[['no_indiv']] <- results.save$no_indiv
  out[['rank_spec']] <- results.save$rank_spec
  
  return(out)
}

# Get paths of all saved files starting with "demo2_"
file_paths <- list.files(path = './results', pattern = "^demo2_0_FALSE_FALSE.*\\.RData$", full.names = TRUE)

# Extract parameters from each file
all_params <- lapply(file_paths, extract_params)

# Convert list of lists to data frame
df <- do.call(rbind, lapply(all_params, data.frame, row.names = NULL))

# Organize the data frame
models <- c('jive', 'slide', 'ajive', 'dcca', 'unifac', 'proposed', 'proposed_subsampling')
cols <- sapply(models, function(x) paste0(x,'_avgF1'))
names(cols) <- NULL
df <- df[, c(cols,
             "phi_max", 
             "SNR1", 
             "SNR2", 
             "no_joint", 
             "no_indiv",
             "rank_spec")]

# Print the data frame
df_short <- df[, cols]
df_short <- cbind(expand.grid(c(90, 60, 45, 30), c(8, 4, 2, 1)), df_short)
# round to 3 decimal places
df_short <- round(df_short, 3)
# kable(df_short[, c('Var2', 'Var1', 'f1_jive', 'f1_slide', 'f1_ajive', 'f1_dcca', 'f1_unifac', 'f1_proposed', 'f1_proposed_subsampling')], 'latex')
