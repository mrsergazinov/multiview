convert_to_latex <- function(df) {
  # Initialize an empty vector to store the LaTeX rows
  latex_rows <- vector("character", ncol(df))
  
  # Loop through each column to create the LaTeX formatted string
  for (i in seq_along(df)) {
    column_name <- colnames(df)[i]
    column_values <- df[[i]]
    
    # Create a LaTeX formatted row with column name and its values
    latex_rows[i] <- paste(column_name, "&", paste(column_values, collapse = " & "), "\\\\")
  }
  
  # Combine all rows into a single LaTeX table format
  latex_table <- paste(latex_rows, collapse = "\n")
  
  return(latex_table)
}