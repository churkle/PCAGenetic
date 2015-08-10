perform_PCA <- function(file, colStart = 1, colEnd = max(count.fields(file, skip = colStart)), 
                        alleleCol = "alleles", numPCs = 3)
{
    pc1ThetaMin <- 30
    pc1ThetaMax <- 150
    pc1PhiMin <- -30
    pc1PhiMax <- 90
    pc2PhiMin <- 90
    pc2PhiMax <- 210
    pc3Max <- 30
  
    mydata <- read.table(file = file, header = TRUE, colClasses = "character")
    bases <- mydata[ colStart:colEnd ]
    pcSorted <- data.frame("PC1" = character(0), "PC2" = character(0),
                           "PC3" = character(0), stringsAsFactors = FALSE)


    for(i in 1:nrow(mydata))
    {
        letterSplit <- strsplit(mydata[ i,alleleCol ], split = "/")
        letterSplit <- sort(letterSplit[[1]])
        baseVect <- character(3)
        baseVect[1] <- paste(letterSplit[1], letterSplit[1], sep = "")
        baseVect[2] <- paste(letterSplit[1], letterSplit[2], sep = "")
        baseVect[3] <- paste(letterSplit[2], letterSplit[2], sep = "")
        
        for(j in colStart:ncol(mydata))
        {
          bases[ i,j - colStart + 1 ] <- match(as.character(mydata[ i,j ]), baseVect) - 1
        }
    }
    
    bases <- data.matrix(bases)
    bases <- t(bases)
    myPCA <- prcomp(bases)
    capture.output( print(myPCA$rotation, print.gap = 1, row.names = TRUE, right = F), file = "loading.txt" )
    capture.output( print(myPCA$x, print.gap = 1, row.names = FALSE, right = F), file = "scores.txt" )
    
    loading <- myPCA$rotation[ ,1:3]
    score <- myPCA$x
    #pcaRegion <- matrix(nrow = nrow(loading), ncol = 3, dimnames = list(NULL, c("PC1", "PC2", "PC3")))
    
    for(k in seq(0, 180, 30))
    {
      bAngle <- (k / 180) * pi
      
      for(p in seq(0, 180, 30))
      {
        cAngle <- (p / 180) * pi
        
        for(g in seq(0, 180, 30))
        {
          numPC1 <- 1
          numPC2 <- 1
          numPC3 <- 1
          
          pc1Col <- character()
          pc2Col <- character()
          pc3Col <- character()
          
          dAngle <- (p / 180) * pi
          bMatrix <- matrix(c(1, 0, 0, 0, cos(bAngle), -sin(bAngle), 0, sin(bAngle), cos(bAngle)), 
                            nrow = 3, ncol = 3)
          cMatrix <- matrix(c(cos(cAngle), 0, sin(cAngle), 0, 1, 0, -sin(cAngle),
                              0, cos(cAngle)), nrow = 3, ncol = 3)
          dMatrix <- matrix(c(cos(dAngle), -sin(dAngle), 0, sin(dAngle), cos(dAngle),
                              0, 0, 0, 1), nrow = 3, ncol = 3)
          dcMatrix <- dMatrix %*% cMatrix
          dcbMatrix <- dcMatrix %*% bMatrix

          newLoading <- loading %*% dcbMatrix
          row_sub <- apply(newLoading, 1, function(row) all(row !=0 ))
          newLoading <- newLoading[row_sub,]
          newMyData <- mydata[row_sub,]
          spherical <- matrix(nrow = nrow(newLoading), ncol = 3, 
                              dimnames = list(NULL ,c("r", "theta", "phi")))
          

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
          
          for(t in 1:nrow(newLoading))
          {
            if(spherical[ t,2 ] >= pc1ThetaMin & spherical[ t,2 ] <= pc1ThetaMax
               & spherical[ t,3 ] >= pc1PhiMin & spherical[ t,3 ] <= pc1PhiMax)
            {
              curNumPC1 <- numPC1
              for(m in 0:(colEnd - colStart))
              {
                pcSorted[ m + curNumPC1,"PC1" ] <- newMyData[ t,colStart + m ]
                numPC1 <- numPC1 + 1
#                 numPC1 <- numPC1 + 1
#                 pc1Col[numPC1] <- mydata[ t,colStart + m - 1 ]
              }
            }
            else if(spherical[ t,2 ] >= pc1ThetaMin & spherical[ t,2 ] <= pc1ThetaMax
                    & spherical[ t,3 ] >= pc2PhiMin & spherical[ t,3 ] <= pc2PhiMax)
            {
              curNumPC2 <- numPC2
              for(n in 0:(colEnd - colStart))
              {
                pcSorted[ n + curNumPC2,"PC2" ] <- newMyData[ t,colStart + n ]
                numPC2 <- numPC2 + 1
#                 numPC2 <- numPC2 + 1
#                 pc1Col[numPC2] <- mydata[ t,colStart + m - 1 ]
              }
            }
            else if(spherical[ t,2 ] <= pc3Max)
            {
              curNumPC3 <- numPC3
              for(o in 0:(colEnd - colStart))
              {
                pcSorted[ o + curNumPC3,"PC3" ] <- newMyData[ t,colStart + o ]
                numPC3 <- numPC3 + 1
#                 numPC3 <- numPC3 + 1
#                 pc3Col[numPC3] <- mydata[ t,colStart + m - 1 ]
              }
            }
          }
        }
      }
    }
    generate_BCS_data(pcSorted)
}

