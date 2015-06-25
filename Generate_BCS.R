generate_BCS <- function(file, colStart = 1, colEnd = max(count.fields(file, skip = colStart)))
{
  
  mydata <- read.table(file = file, header = TRUE)[, colStart:colEnd]
  bcsData <- data.frame("Strain" = character(0), stringsAsFactors = FALSE)
  
  for(i in 1:ncol(mydata))
  {
    bcsData[i, "Strain"] <- colnames(mydata)[i]
    baseCount <- new.env(hash = TRUE)
    total <- 0
    
    for(j in 1:ncol(mydata))
    {
      
      if(exists(as.character(mydata[ i,j ]), where = baseCount))
      {
        baseCount[[ as.character(mydata[ i,j ]) ]] <- 
          baseCount[[ as.character(mydata[ i,j ]) ]] + 1
      }
      else
      {
        if(!(as.character(mydata[ i,j ]) %in% colnames(bcsData)))
        {
          bcsData[[ as.character(mydata[ i,j ]) ]] <- 0
        }
        baseCount[[ as.character(mydata[ i,j ]) ]] <- 1
      }
      total <- total + 1
    }
    
    for(k in ls(baseCount))
    {
      bcsData[ i,k ] <- baseCount[[ as.character(k) ]]
    }
  }
  bcsData
}