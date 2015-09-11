#'Perform PCA
#'
#'This function performs PCA on a hmp file and 
#'compares the PC BCS data to the PCA scores
#'@param file The hmp file that contains genomic data
#'@param colStart The column to start reading data from
#'@param colEnd The last column to read data from
#'@param alleleCol The name of the column that contains the allele for the row
#'@param numPCs The number of PCs to generate
#'@examples
#'perform_PCA("test_files/H1k_0_01.hmp", 12)

perform_PCA <- function(file, colStart = 1, colEnd = max(count.fields(file, skip = colStart)), 
                        alleleCol = "alleles", numPCs = 3)
{
    #Angle thresholds for each PC
    pc1ThetaMin <- 30
    pc1ThetaMax <- 150
    pc1PhiMin <- -30
    pc1PhiMax <- 90
    pc2PhiMin <- 90
    pc2PhiMax <- 210
    pc3Max <- 30
  
    #Read data from file and put it into a data frame
    mydata <- read.table(file = file, header = TRUE, colClasses = "character")
    bases <- mydata[ colStart:colEnd ]
    
    #Data frame for holding results
    bcsData <- data.frame("A" = numeric(0), "T" = numeric(0), 
                          "G" = numeric(0), "C" = numeric(0), stringsAsFactors = FALSE)
    
    #Data frame for holding correlations
    correlations <- data.frame("Rotation" = character(0), "R1" = numeric(0),
                               "R2" = numeric(0), "R3" = numeric(0),
                               "R Total" = numeric(0), stringsAsFactors = FALSE)
    rotationNum <- 1

    #Loop through every row in the data
    for(i in 1:nrow(mydata))
    {
        #Split multiple letter data into separate characters and sort them 
        #alphabetically
        letterSplit <- strsplit(mydata[ i,alleleCol ], split = "/")
        letterSplit <- sort(letterSplit[[1]])
        
        #Find the 3 variations of the possible combinations and arrange them
        #alphabetically
        baseVect <- character(3)
        baseVect[1] <- paste(letterSplit[1], letterSplit[1], sep = "")
        baseVect[2] <- paste(letterSplit[1], letterSplit[2], sep = "")
        baseVect[3] <- paste(letterSplit[2], letterSplit[2], sep = "")
        
        #Iterate through the data columns and convert the character combination to a number
        for(j in colStart:ncol(mydata))
        {
          bases[ i,j - colStart + 1 ] <- match(as.character(mydata[ i,j ]), baseVect) - 1
        }
    }
    
    #Transpose data and perform PCA on it
    bases <- data.matrix(bases)
    bases <- t(bases)
    myPCA <- prcomp(bases)
    
    loading <- myPCA$rotation[ ,1:3]
    score <- myPCA$x
    
    #Generate Euler Rotations
    for(k in seq(0, 180, 30))
    {
      bAngle <- (k / 180) * pi
      
      for(p in seq(0, 180, 30))
      {
        cAngle <- (p / 180) * pi
        
        for(g in seq(0, 180, 30))
        {
          #Set all BCS totals to 0
          pc1Total <- 0
          pc2Total <- 0
          pc3Total <- 0
          pc1ATotal <- 0
          pc1TTotal <- 0
          pc1GTotal <- 0
          pc1CTotal <- 0
          pc2ATotal <- 0
          pc2TTotal <- 0
          pc2GTotal <- 0
          pc2CTotal <- 0
          pc3ATotal <- 0
          pc3TTotal <- 0
          pc3GTotal <- 0
          pc3CTotal <- 0
          
          #Rotate the matrix by the Euler rotations
          dAngle <- (p / 180) * pi
          bMatrix <- matrix(c(1, 0, 0, 0, cos(bAngle), -sin(bAngle), 0, sin(bAngle), cos(bAngle)), 
                            nrow = 3, ncol = 3)
          cMatrix <- matrix(c(cos(cAngle), 0, sin(cAngle), 0, 1, 0, -sin(cAngle),
                              0, cos(cAngle)), nrow = 3, ncol = 3)
          dMatrix <- matrix(c(cos(dAngle), -sin(dAngle), 0, sin(dAngle), cos(dAngle),
                              0, 0, 0, 1), nrow = 3, ncol = 3)
          dcMatrix <- dMatrix %*% cMatrix
          dcbMatrix <- dcMatrix %*% bMatrix

          #Find the new loading of the matrix and remove zero values
          newLoading <- loading %*% dcbMatrix
          row_sub <- apply(newLoading, 1, function(row) all(row !=0 ))
          newLoading <- newLoading[row_sub,]
          newMyData <- mydata[row_sub,]
          
          #Create matrix to hold spherical coordinates of rotated values
          spherical <- matrix(nrow = nrow(newLoading), ncol = 3, 
                              dimnames = list(NULL ,c("r", "theta", "phi")))
          
          #Calculate spherical coordinates of the newly rotated values
          for(l in 1:nrow(newLoading))
          {
            x <- newLoading[ l,1 ]
            y <- newLoading[ l,2 ]
            z <- newLoading[ l,3 ]
            
            r <- sqrt(x^2 + y^2 + z^2)
            theta <- acos(z/r)
            phi <- atan(y/x)
            
            spherical[ l,1 ] <- r
            spherical[ l,2 ] <- theta
            spherical[ l,3 ] <- phi
          }
          
          #Find BCS for each individual PC
          #REWRITE THIS PART
          for(t in 1:nrow(newLoading))
          {
            #PC1
            if(spherical[ t,2 ] >= pc1ThetaMin & spherical[ t,2 ] <= pc1ThetaMax
               & spherical[ t,3 ] >= pc1PhiMin & spherical[ t,3 ] <= pc1PhiMax)
            {
                for(q in colStart:colEnd)
                {
                    letters <- strsplit(as.character(newMyData[ t,q ]), split = "")
                    for(g in 1:length(letters[[1]]))
                    {
                      if(letters[[1]][g] == "A")
                      {
                        pc1ATotal <- pc1ATotal + 1
                      }
                      else if(letters[[1]][g] == "T")
                      {
                        pc1TTotal <- pc1TTotal + 1
                      }
                      else if(letters[[1]][g] == "G")
                      {
                        pc1GTotal <- pc1GTotal + 1
                      }
                      else if(letters[[1]][g] == "C")
                      {
                        pc1CTotal <- pc1CTotal + 1
                      }
                      pc1Total <- pc1Total + 1
                    }
                }
            }
            #PC2
            else if(spherical[ t,2 ] >= pc1ThetaMin & spherical[ t,2 ] <= pc1ThetaMax
                    & spherical[ t,3 ] >= pc2PhiMin & spherical[ t,3 ] <= pc2PhiMax)
            {
              for(q in colStart:colEnd)
              {
                letters <- strsplit(as.character(newMyData[ t,q ]), split = "")
                for(g in 1:length(letters[[1]]))
                {
                  if(letters[[1]][g] == "A")
                  {
                    pc2ATotal <- pc2ATotal + 1
                  }
                  else if(letters[[1]][g] == "T")
                  {
                    pc2TTotal <- pc2TTotal + 1
                  }
                  else if(letters[[1]][g] == "G")
                  {
                    pc2GTotal <- pc2GTotal + 1
                  }
                  else if(letters[[1]][g] == "C")
                  {
                    pc2CTotal <- pc2CTotal + 1
                  }
                  pc2Total <- pc2Total + 1
                }
              }
            }
            #PC3
            else if(spherical[ t,2 ] <= pc3Max)
            {
              for(q in colStart:colEnd)
              {
                letters <- strsplit(as.character(newMyData[ t,q ]), split = "")
                for(g in 1:length(letters[[1]]))
                {
                  if(letters[[1]][g] == "A")
                  {
                    pc3ATotal <- pc3ATotal + 1
                  }
                  else if(letters[[1]][g] == "T")
                  {
                    pc3TTotal <- pc3TTotal + 1
                  }
                  else if(letters[[1]][g] == "G")
                  {
                    pc3GTotal <- pc3GTotal + 1
                  }
                  else if(letters[[1]][g] == "C")
                  {
                    pc3CTotal <- pc3CTotal + 1
                  }
                  pc3Total <- pc3Total + 1
                }
              }
            }
          }
          
          #Fill data frame with BCS results from individual PCs
          bcsData[[ "PC1","A" ]] <- pc1ATotal / pc1Total
          bcsData[[ "PC1","T" ]] <- pc1TTotal / pc1Total
          bcsData[[ "PC1","G" ]] <- pc1GTotal / pc1Total
          bcsData[[ "PC1","C" ]] <- pc1CTotal / pc1Total
          bcsData[ "PC1","Total" ] <- pc1Total
          bcsData[[ "PC2","A" ]] <- pc2ATotal / pc2Total
          bcsData[[ "PC2","T" ]] <- pc2TTotal / pc2Total
          bcsData[[ "PC2","G" ]] <- pc2GTotal / pc2Total
          bcsData[[ "PC2","C" ]] <- pc2CTotal / pc2Total
          bcsData[ "PC2","Total" ] <- pc2Total
          bcsData[[ "PC3","A" ]] <- pc3ATotal / pc3Total
          bcsData[[ "PC3","T" ]] <- pc3TTotal / pc3Total
          bcsData[[ "PC3","G" ]] <- pc3GTotal / pc3Total
          bcsData[[ "PC3","C" ]] <- pc3CTotal / pc3Total
          bcsData[ "PC3","Total" ] <- pc3Total
          
          #This needs to be fixed
          r1 <- cor(x = bcsData[["PC1", "A"]], y = myPCA$x[["PC1"]])
          r2 <- cor(x = bcsData[["PC2", "A"]], y = myPCA$x[["PC2"]])
          r3 <- cor(x = bcsData[["PC3", "A"]], y = myPCA$x[["PC3"]])
          correlations[[""]]
          rotationNum <- rotationNum + 1
        }
      }
    }
    write(bcsData, file = "PCBCData.txt")
    bcsData
}

