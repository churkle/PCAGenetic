#' Generate BCS Function
#' 
#' This function generates BCS for an hmp file and generates a histogram
#' @param file The hmp file to generate BCS for
#' @param colStart The column to start counting base data from
#' @param colEnd The last column to count
#' @examples
#' generate_BCS("test_files/H1k_0_01.hmp", 12)

generate_BCS <- function(file, colStart = 1, colEnd = max(count.fields(file, skip = colStart)))
{
  #Read data from file and put it into a data frame
  mydata <- read.table(file = file, header = TRUE)[, colStart:colEnd]
  
  #Create data frame to put results into
  bcsData <- data.frame("Strain" = character(0), "A" = numeric(0), "T" = numeric(0), 
                        "G" = numeric(0), "C" = numeric(0), stringsAsFactors = FALSE)

  #Go through every column in the data
  for(i in 1:ncol(mydata))
  {
    #Give the result data frame the correct column names
    bcsData[i, "Strain"] <- colnames(mydata)[i]
    bcsData[i, !names(bcsData) %in% "Strain"] <- 0
    
    #The total of each base
    total <- 0
    aCount <- 0
    tCount <- 0
    gCount <- 0
    cCount <- 0
    
    #Go through each row in the data
    for(j in 1:nrow(mydata))
    {
      #Split double-based data into separate characters
      letters <- strsplit(as.character(mydata[ j,i ]), split = "")
      
      #Loop through each character in the double-based data
      for(g in 1:length(letters[[1]]))
      {
          #Increment the total for each base character if it matches
          if(letters[[1]][g] == "A")
          {
              aCount <- aCount + 1
          }
          else if(letters[[1]][g] == "T")
          {
              tCount <- tCount + 1
          }
          else if(letters[[1]][g] == "G")
          {
              gCount <- gCount + 1
          }
          else if(letters[[1]][g] == "C")
          {
              cCount <- cCount + 1
          }
          
          #Increment the total amount of characters
          total <- total + 1
      }
    
    }
        
    #Find the BCS by dividing base frequency by the total amount
    bcsData[[ i,"A" ]] <- aCount / total
    bcsData[[ i,"T" ]] <- tCount / total
    bcsData[[ i,"G" ]] <- gCount / total
    bcsData[[ i,"C" ]] <- cCount / total
    
    bcsData[ i,"Total" ] <- total
    
  }
  
  #Write the output data to a text file
  capture.output( print(bcsData, print.gap = 1, row.names = FALSE, right = F), 
                  file = "output.txt" )
  
  #Generate a histogram of the BCS data
  hist(bcsData[ ,"A" ], xlab = "A", breaks = 200)
}