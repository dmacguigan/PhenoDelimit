# read in CLUMPP results for permuted runs from miscfiles
read_clumpp_results_permuted <- function(wd, perc.var, model.numbers, permutations){
  setwd(wd)
  c1=numeric()
  c2=numeric()
  c3=numeric()
  c4=numeric()
  for(i in 1:length(perc.var)){
    for(j in 1:length(model.numbers)){
      for(p in 1:permutations){
        fileName <- Sys.glob(paste("m", model.numbers[j], "_perVar-", perc.var[i], ".p", p, ".miscfile", sep=""))
        temp <- grep("The highest value of H' is:", readLines(fileName), value=TRUE)
        H <- as.numeric(strsplit(temp, " ")[[1]][7]) # get only H' number
        c1 <- c(c1, model.numbers[j])
        c2 <- c(c2, perc.var[i])
        c3 <- c(c3, H)
        c4 <- c(c4, p)
      }
    }
  }
  
  # create data frame with results
  H_df <- data.frame(model=c1, perc.pca=c2, H=c3, permutation=c4)
  
  return(H_df)
}

