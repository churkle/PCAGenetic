perform_PCA <- function(file, colStart = 1, colEnd = max(count.fields(file, skip = colStart)), alleleCol = "alleles")
{
    mydata <- read.table(file = file, header = TRUE, colClasses = "character")
    bases <- mydata[ colStart:colEnd ]

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
    prcomp(bases)
}