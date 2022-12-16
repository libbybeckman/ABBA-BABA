Pattern_By_Site <- function(alignment = "alignment.fasta"){
  
  # Purpose of function is to characterize each site in a fasta file as "ABBA","BABA" or "Uninformative"
  # Requires "ape" and "seqinr" libraries to be loaded
  # Run source code before using in command line
  # To run in the console:
  # Result <- Pattern_By_Site(alignment="/path to alignment.fasta")
  # Result is a data.frame
  
  
  # Load in empirical data in fasta format
  alignment <- read.alignment(alignment, format = "fasta")                         #  read in the alignment
  alignment.matrix <- matrix(, length(alignment$nam), nchar(alignment$seq[[1]]))    #  make a matrix for the alignment
  for(i in 1:length(alignment$nam)){
    alignment.matrix[i, ] <- unlist(strsplit(alignment$seq[[i]], ""))               #  fill in the matrix
  }
  
  abba <- 0                                                                         #  set up my variables
  baba <- 0                                                                         #  set up my variables
  position <- 0
  Sites <- ncol(alignment.matrix)
  df <- data.frame(matrix(ncol = 2, nrow = Sites))
  colnames(df) <- c("Site","Pattern")
  for (i in 1:ncol(alignment.matrix)){
    position <- position + 1
    df[i,1] <- position
    if(length(unique(alignment.matrix[, i])) == 2){                                 #  unique(c(p1,p2,p3,o))==2 aka biallelic
      if(alignment.matrix[1, i] != alignment.matrix[2, i]){                         #  p1 != p2   aka different resolutions in p1 and p2
        if(alignment.matrix[4, i] != alignment.matrix[3, i]){                       #  o != p3    durand says "less likely pattern due to seq. errors
          if(alignment.matrix[3, i] == alignment.matrix[1, i]) {
            baba <- baba + 1                                          #  add to the count of baba sites
            df[i,2] <- "BABA"}                                        # Add BABA to pattern column
          if(alignment.matrix[2, i] == alignment.matrix[3, i]) {
            abba <- abba + 1                                          #  add to the count of abba sites
            df[i,2] <- "ABBA"}                                        # Add ABBA to pattern column
              } 
       }
    }
    if (is.na(df[i,2])) {
      df[i,2] <- "Uninformative"
    }    
  }
  return(df)
  }