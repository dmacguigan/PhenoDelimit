# plot discriminant axis variable loadings and write table
discriminant_loading <- function(wd, clumpp.wd,
                                 best.perc.var, best.model.number,
                                 plot.type, plot.width, plot.height,
                                 axis){

  # read DAPC variable contributions from file
  var.contr <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_var.contr.txt", sep=""),
                          header = TRUE)

  # read DAPC variable loadings from file
  var.load <- read.table(file = paste(clumpp.wd, "/m", best.model.number, "_perVar-", best.perc.var, "_var.load.txt", sep=""),
                         header=TRUE)

  setwd(wd)

  if(plot.type =="svg"){
    svglite(file = paste("DA_", axis, "_loading.svg", sep=""), width = plot.width, height = plot.height)
    loadingplot(var.contr, axis=axis)
    dev.off()
  } else if(plot.type =="pdf"){
    pdf(file = paste("DA_", axis, "_loading.pdf", sep=""), width = plot.width, height = plot.height)
    loadingplot(var.contr, axis=axis)
    dev.off()
  }

  d <- var.contr[,axis]
  d2 <- var.load[,axis]
  t <- cbind(d,d2)
  t <- as.data.frame(t)
  row.names(t) <- row.names(var.contr)
  t <- t[order(t$d, decreasing = TRUE),]
  colnames(t) <- c("contribution", "loading")
  print(head(t))
  print(str(t))
  write.table(x = t, file = paste("DA_", axis, "_loading.txt", sep=""), quote=FALSE, row.names = TRUE)
}


