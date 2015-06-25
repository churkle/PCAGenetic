generate_BCS <- function(file, colStart = 1, colEnd = max(count.fields(file, skip = colStart)))
{
  
  mydata <- read.table(file = file, header = TRUE)[, colStart:colEnd]
  bcsData <- data.frame("Strain" = character(0), stringsAsFactors = FALSE)
  
  for(i in 1:ncol(mydata))
  {
    bcsData[i, "Strain"] <- colnames(mydata)[i]
    bcsData[i, !names(bcsData) %in% "Strain"] <- 0
    baseCount <- new.env(hash = TRUE)
    total <- 0
    
    for(j in 1:nrow(mydata))
    {
      letters <- strsplit(as.character(mydata[ j,i ]), split = "")
      
      for(g in i:length(letters))
      {
          if(exists(letters[g], where = baseCount))
          {
            baseCount[[ as.character(letters[g]) ]] <- 
              baseCount[[ as.character(letters[g]) ]] + 1
          }
          else
          {
            if(!(letters[g]) %in% colnames(bcsData))
            {
              bcsData[[ letters[g] ]] <- 0
            }
            baseCount[[ letters[g] ]] <- 1
          }
          total <- total + 1
      }
      
    }
    
    for(k in ls(baseCount))
    {
      bcsData[ i,k ] <- baseCount[[ as.character(k) ]] / total
    }
    
    bcsData[ i,"Total" ] <- total
    
  }
  
  bcsData
}