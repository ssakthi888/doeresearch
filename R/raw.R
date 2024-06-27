rawdf=function(input_filename,output_filename){
  #load libraries
  library(reshape2)
  library(tidyr)
  library(dplyr)
  
  # Read the CSV file
  data <- read.csv(input_filename,check.names = FALSE)
  
  # Get the number of columns in the original data frame
  num_columns <- ncol(data)
  
  # Initialize a counter for the dataframes
  df_count <- 0
  
  # Create data frames with four columns each and assign specific names
  for (i in seq(1, num_columns, by = 4)) {
    start_col <- i
    end_col <- min(i + 3, num_columns)
    
    # Create a new data frame with the selected columns
    df <- data[, start_col:end_col]
    
    # Increment the dataframe counter
    df_count <- df_count + 1
    # Assign a specific name to the data frame
    assign(paste("df", df_count, sep = ""), df)
  }
  
  # Print the total number of dataframes generated
  cat("Total number of dataframes:", df_count, "\n")
  
  # Initialize an empty list to store the melted dataframes
  melted_dfs <- list()
  
  # Iterate over the list of dataframe names
  for (i in 1:df_count) {
    df_name <- paste("df", i, sep = "")
    first_col_title <- colnames(get(df_name))[1]
    
    # Use melt function from reshape2, retaining the first column
    df_f <- melt(get(df_name), id.vars = first_col_title, variable.name = "Rep", value.name = paste(first_col_title, "_V", sep = ""))
    
    # Separate the first column into "Main" and "Sub"
    df_f <- separate(df_f, first_col_title, into = c("Main", "Sub"), sep = 2)
    
    # Remove "_V" from the variable names
    colnames(df_f) <- sub("_V$", "", colnames(df_f))
    
    
    # Assign the melted dataframe to a list
    melted_dfs[[df_name]] <- df_f
  }
  
  
  # Save the combined dataframe to a CSV file
  a=write.csv(melted_dfs, file = "1245.csv", row.names = FALSE)
  a=read.csv("1245.csv",check.names = FALSE)
  file.remove("1245.csv")
  
  # Assuming your dataframe is named 'a'
  old_colnames <- colnames(a)
  
  # Remove the prefixes (df1., df2., etc.)
  new_colnames <- sub("^df[0-9]+\\.", "", old_colnames)
  
  # Set the new column names
  colnames(a) <- new_colnames
  
  # Remove duplicate columns
  a <- a[, !duplicated(colnames(a))]
  
  # Save the combined dataframe to a CSV file
  a=write.csv(a, file = output_filename, row.names = FALSE)
  
  print("code finished")
}
