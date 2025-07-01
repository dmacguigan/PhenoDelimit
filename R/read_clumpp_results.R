#' Summarize CLUMPP results
#' 
#' Read in results from CLUMPP miscfiles, return data frame of H' values for each delimitation model
#' 
#' @param wd working directory to store CLUMPP files, should contain a CLUMPP_permuted directory created by the dapc_clumpp_permuted function
#' @param model.numbers vector containing delimitation model numbers
#' @param perc.var vector containing cumulative percentages of variance to retain for discriminant analyses
#' 
#' @export
#' 
read_clumpp_results <- function(wd, perc.var, model.numbers, apriori=FALSE, clust.method="kmeans"){
  # input error handling
	if(any(perc.var <= 0) || (any(perc.var > 100))){
		stop("perc.var values must be greater than 0 and less than 100")
	}

  if(clust.method == "kmeans" && apriori == FALSE){
	  setwd(paste0(wd, "/CLUMPP"))
	  c1=numeric()
	  c2=numeric()
	  c3=numeric()
	  for(i in 1:length(perc.var)){
		for(j in 1:length(model.numbers)){
		  fileName <- Sys.glob(paste("m", model.numbers[j], "_perVar-", perc.var[i], ".miscfile", sep=""))
		  temp <- grep("The highest value of H' is:", readLines(fileName), value=TRUE)
		  H <- as.numeric(strsplit(temp, " ")[[1]][7]) # get only H' number
		  c1 <- c(c1, model.numbers[j])
		  c2 <- c(c2, perc.var[i])
		  c3 <- c(c3, H)
		}
	  }

	  # create data frame with results
	  H_df <- data.frame(model=c1, perc.pca=c2, H=c3)

	  # which model is the max for each plotting group?
	  maxes <- aggregate(H_df$H, by = list(H_df$perc.pca), max)
	  H_df$max <- ifelse((H_df$H %in% maxes$x), TRUE, FALSE)

	  return(H_df)
  }else if(clust.method == "randomforest" && apriori == FALSE){
	  setwd(paste0(wd, "/CLUMPP"))
	  c1=numeric()
	  c2=numeric()
	  c3=numeric()
	  for(i in 1:length(perc.var)){
		for(j in 1:length(model.numbers)){
		  fileName <- Sys.glob(paste("m", model.numbers[j], "_perVar-", perc.var[i], ".RF.miscfile", sep=""))
		  temp <- grep("The highest value of H' is:", readLines(fileName), value=TRUE)
		  H <- as.numeric(strsplit(temp, " ")[[1]][7]) # get only H' number
		  c1 <- c(c1, model.numbers[j])
		  c2 <- c(c2, perc.var[i])
		  c3 <- c(c3, H)
		}
	  }

	  # create data frame with results
	  H_df <- data.frame(model=c1, perc.pca=c2, H=c3)

	  # which model is the max for each plotting group?
	  maxes <- aggregate(H_df$H, by = list(H_df$perc.pca), max)
	  H_df$max <- ifelse((H_df$H %in% maxes$x), TRUE, FALSE)

	  return(H_df)
  }else if(clust.method == "randomforest" && apriori == TRUE){
	  setwd(paste0(wd, "/CLUMPP"))
	  c1=numeric()
	  c2=numeric()
	  c3=numeric()
	  for(i in 1:length(perc.var)){
		for(j in 1:length(model.numbers)){
		  fileName <- Sys.glob(paste("m", model.numbers[j], "_perVar-", perc.var[i], ".RF.supervised.miscfile", sep=""))
		  temp <- grep("The highest value of H' is:", readLines(fileName), value=TRUE)
		  H <- as.numeric(strsplit(temp, " ")[[1]][7]) # get only H' number
		  c1 <- c(c1, model.numbers[j])
		  c2 <- c(c2, perc.var[i])
		  c3 <- c(c3, H)
		}
	  }

	  # create data frame with results
	  H_df <- data.frame(model=c1, perc.pca=c2, H=c3)

	  # which model is the max for each plotting group?
	  maxes <- aggregate(H_df$H, by = list(H_df$perc.pca), max)
	  H_df$max <- ifelse((H_df$H %in% maxes$x), TRUE, FALSE)

	  return(H_df)
  }
}


