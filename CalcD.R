CalcD <- function(alignment = "alignment.fasta", 
                  sig.test="N",                                                    # options are "N", "B", "J"
                  block.size = 1000,                                                # size of blocks to drop in jacknife
                  replicate=1000){
  # this function is used regardless of approach
  d.calc <- function(alignment){
    abba <- 0                                                                         #  set up my variables
    baba <- 0                                                                         #  set up my variables
    for(i in 1:ncol(alignment)){                                               #  run through all sites
      if(length(unique(alignment[, i])) == 2){                                 #  unique(c(p1,p2,p3,o))==2 aka biallelic
        if(alignment[1, i] != alignment[2, i]){                         #  p1 != p2   aka different resolutions in p1 and p2
          if(alignment[4, i] != alignment[3, i]){                       #  o != p3    durand says "less likely pattern due to seq. errors
            if(alignment[3, i] == alignment[1, i]) {baba <- baba + 1}   #  add to the count of baba sites
            if(alignment[2, i] == alignment[3, i]) {abba <- abba + 1}   #  add to the count of abba sites
          } 
        }
      }
    }
    d <- (abba - baba) / (abba + baba)   #what its all about
    results <- list()
    results[[1]] <- d
    results[[2]] <- abba
    results[[3]] <- baba
    return(results)
  }
  
  #### Test of empirical data
  alignment <- read.alignment(alignment, format = "fasta")                         #  read in the alignment
  alignment.matrix <- matrix(, length(alignment$nam), nchar(alignment$seq[[1]]))    #  make a matrix for the alignment
  for(i in 1:length(alignment$nam)){
    alignment.matrix[i, ] <- unlist(strsplit(alignment$seq[[i]], ""))               #  fill in the matrix
  }
  results <- d.calc(alignment.matrix)
  D <- results[[1]]
  if(is.nan(D)) D <- 0
  abba <- results[[2]]
  baba <- results[[3]]
  
  
  
  ## THIS SECTION WILL CALCULATE THE P-VAL BASED ON BOOTSTRAPPING
  ## SITES ARE SAMPLED WITH REPLACEMENT TO MAKE A NEW DATASET OF
  ## OF EQUAL SIZE TO THE ORIGINAL DATASET THIS ALLOWS US TO CALCULATE
  ## THE STANDARD DEVIATION AND THUS A Z SCORE.
  if(sig.test=="B"){
    sim.d<-vector()
    foo <- ncol(alignment.matrix)
    sim.matrix<-matrix(,4,foo)
    cat("\nperforming bootstrap")
    for(k in 1:replicate){
      if(k/(replicate/100) == round(k/(replicate/100))) cat(".")
      sim.matrix[1:4,1:foo] <- alignment.matrix[1:4, sample(1:foo, replace=T)]
      results <- d.calc(sim.matrix)
      sim.d[k] <- results[[1]]
    }
    sim.d[is.nan(sim.d)] <- 0
    z <- abs(D/sd(sim.d))
    sd.D <- sd(sim.d)
    new.pval <- 2 * (1 - pnorm(z))
    ## NOW WE MAKE THE OUTPUTS  
    #d <- results[[1]]
    #print(d)
    Output <- vector()
    Output[1] <- D            			# d
    Output[2] <- abba       				# abba
    Output[3] <- baba       				# baba
    Output[4] <- z                  # z score
    Output[5] <- sd.D
    Output[6] <- new.pval           # p-value
    Output[7] <- foo-1
    
    #cat("\nSites in alignment =", ncol(alignment.matrix))
    #cat("\nNumber of sites with ABBA pattern =", abba)
    #cat("\nNumber of sites with BABA pattern =", baba)
    #cat("\n\nD raw statistic / Z-score = ", d, " / ", z)
    #cat("\n\nResults from ", replicate, "bootstraps")
    #cat("\nSD D statistic =", sd(sim.d))
    #cat("\nP-value (that D=0) = ",new.pval) #after Eaton and Ree 2013 
  }
  
  ## THIS SECTION WILL CALCULATE THE P-VAL BASED ON JACKKNIFING
  ## THE DATA IS REANALYZED WHILE DROPPING A PORTION OF THE DATA
  ## THE SIZE OF THE DROPPED PORTION IS DETERMINED BY THE BLOCK SIZE
  ## ARGUMENT THIS PROCEDURE ALLOWS US TO CALCULATE
  ## THE STANDARD DEVIATION AND THUS A Z SCORE. THIS APPROACH IS PARTICULARLY 
  ## IMPORTANT WHEN WE WILL BE USING DATA WHERE THE SNPs MAY BE IN LINKAGE WITH
  ## ONE ANOTHER
  if(sig.test=="J"){
    
    #first lets test whether we are limited in the number of reps by block size
    max.rep <- ncol(alignment.matrix) - block.size
    if(block.size >= (ncol(alignment.matrix)/2)){
      stop(call. = F, paste("\nWith a block size of", block.size, 
                            "and an alignment of", ncol(alignment.matrix), 
                            "sites \nsome sites would never be included in the \nanalysis",
                            "\n\nThe maximum block size is 1/2 the alignment length"))
    }
    if(max.rep < replicate){
      stop(call. = F, paste("\nWith a block size of", block.size, 
                            "and an alignment of", ncol(alignment.matrix), 
                            "sites", replicate, "replicates\nare not possible"))
    }
    
    
    if(max.rep >= replicate){
      drop.pos <- seq.int(from=1, to=(max.rep-1), length.out=replicate)
      replicate2 <- replicate
    }
    sim.d<-vector()
    foo <- ncol(alignment.matrix)
    sim.matrix<-matrix(,4,foo-block.size)
    cat("\nperforming jackknife")
    for(k in 1:replicate2){  
      if(k/2 == round(k/2)) cat(".")
      sim.matrix[1:4,1:(foo-block.size-1)] <-alignment.matrix[1:4, -drop.pos[k]:-(drop.pos[k]+block.size)]
      results <- d.calc(sim.matrix)
      sim.d[k] <- results[[1]]
    }
    sim.d[is.nan(sim.d)] <- 0
    z <- abs(d/sd(sim.d))
    new.pval <- 2 * (1 - pnorm(z))
    ## NOW WE MAKE THE OUTPUTS  
    #cat("\nSites in alignment =", ncol(alignment.matrix))
    #cat("\nNumber of sites with ABBA pattern =", abba)
    #cat("\nNumber of sites with BABA pattern =", baba)
    #cat("\nD raw statistic", d)
    #cat("\nZ-score = ", z)
    #cat("\n\nResults from", replicate2, "jackknifes with block size of", block.size)
    #cat("\nSD D statistic =", sd(sim.d))
    #cat("\nP-value (that D=0) = ",new.pval) #after Eaton and Ree 2013 
  }
  if(sig.test=="N"){
    cat("\nSites in alignment =", ncol(alignment.matrix))
    cat("\nNumber of sites with ABBA pattern =", abba)
    cat("\nNumber of sites with BABA pattern =", baba)
    cat("\n\nD raw statistic = ", D)
    Output <- D
  }
  return(Output)
}