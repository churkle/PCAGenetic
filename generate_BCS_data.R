generate_BCS_data <- function(mydata)
{
  bcsData <- data.frame("Strain" = character(0), stringsAsFactors = FALSE)
  
  for(i in 1:ncol(mydata))
  {
    bcsData[i, "Strain"] <- colnames(mydata)[i]
    bcsData[i, !names(bcsData) %in% "Strain"] <- 0
    baseCount <- list()
    total <- 0
    
    for(j in 1:nrow(mydata))
    {
      letters <- strsplit(as.character(mydata[ j,i ]), split = "")
      
      for(g in 1:length(letters[[1]]))
      {
        if(exists(as.character(letters[[1]][g]), where = baseCount) && !is.na(letters[[1]][g]))
        {
          if(!(as.character(letters[[1]][g])) %in% colnames(bcsData))
          {
            bcsData[[ as.character(letters[[1]][g]) ]] <- 0
          }
          
          baseCount[[ as.character(letters[[1]][g]) ]] <- 
            baseCount[[ as.character(letters[[1]][g]) ]] + 1
        }
        else
        {
          if(!(as.character(letters[[1]][g])) %in% colnames(bcsData) && 
             !is.na(letters[[1]][g]))
          {
            bcsData[[ as.character(letters[[1]][g]) ]] <- 0
          }
          
          baseCount[[ as.character(letters[[1]][g]) ]] <- 1
          
        }
        total <- total + 1
      }
      
    }
    
    for(k in names(baseCount))
    {
        bcsData[[ i,k ]] <- baseCount[[ as.character(k)]] / total
    }
    
    bcsData[ i,"Total" ] <- total
    
  }
  
  capture.output( print(bcsData, print.gap = 1, row.names = FALSE, right = F), file = "output.txt" )
  bcsData
}