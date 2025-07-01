#' Summarize permuted CLUMPP results
#' 
#' Read in results from CLUMPP miscfiles, return data frame of H' values for each delimitation model
#' 
#' @param wd working directory to store CLUMPP files, should contain a CLUMPP_permuted directory created by the dapc_clumpp_permuted function
#' @param model.numbers vector containing delimitation model numbers
#' @param perc.var vector containing cumulative percentages of variance to retain for discriminant analyses
#' @param permutations number of permutations performed
#' 
#' @export
#' 
read_clumpp_results_permuted <- function(wd, perc.var, model.numbers, permutations){
  # input error handling
	if(any(perc.var <= 0) || (any(perc.var > 100))){
		stop("perc.var values must be greater than 0 and less than 100")
	}
  if(permutations <=0 ){
		stop("number of permutations must be greater than 0")
	}

  setwd(paste0(wd, "/CLUMPP_permuted"))
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

