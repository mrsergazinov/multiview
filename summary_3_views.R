library(knitr)
source('src/utils.R')

extract_params <- function(file_path) {
  load(file_path)
  
  out <- list()
  models <- names(results.save$results)
  for (model in models) {
    table <- results.save$results[[model]]
    table[['precision.joint']] <- 1 - table$fdr.Pjoint
    table[['precision.indiv1']] <- 1 - table$fdr.Pindiv1
    table[['precision.indiv2']] <- 1 - table$fdr.Pindiv2
    table[['precision.indiv3']] <- 1 - table$fdr.Pindiv3
    table[['avg.precision']] <- (table[['precision.joint']] + 
                                   table[['precision.indiv1']] + 
                                   table[['precision.indiv2']] + 
                                   table[['precision.indiv3']]) / 4
    table[['avg.tpr']] <- (table[['tpr.Pjoint']] + 
                             table[['tpr.Pindiv1']] + 
                             table[['tpr.Pindiv2']] + 
                             table[['tpr.Pindiv3']]) / 4
    table[['f1.avg']] <- 2 * table$avg.precision * table$avg.tpr / (table$avg.precision + table$avg.tpr)
    out[[paste0(model, '_avgF1')]] <- mean(table$f1.avg)
    out[[paste0(model, '_stdF1')]] <- sd(table$f1.avg)
  }
  
  # Extract desired parameters
  out[['phi_max']] <- results.save$phi_max
  out[['SNR1']] <- results.save$SNR1
  out[['SNR2']] <- results.save$SNR2
  out[['SNR3']] <- results.save$SNR3
  out[['no_joint']] <- results.save$no_joint
  out[['no_indiv']] <- results.save$no_indiv
  out[['rank_spec']] <- results.save$rank_spec
  
  return(out)
}

# Get paths of all saved files starting with "demo2_"
file_paths <- list.files(path = './results', pattern = "^demo3_1_FALSE_FALSE.*\\.RData$", full.names = TRUE)

# Extract parameters from each file
all_params <- lapply(file_paths, extract_params)

# Convert list of lists to data frame
df <- do.call(rbind, lapply(all_params, data.frame, row.names = NULL))

# Round to 3 decimal places
df <- round(df, 3)

# Apply the function and print the result
latex_output <- convert_to_latex(df)
cat(latex_output)