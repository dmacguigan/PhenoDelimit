# plot CLUMPP H' values for each delimitaiton model
plot_clumpp_results <- function(wd, clumpp.data, colors, plot.name, plot.type, plot.width, plot.height){
  clumpp.dataBest <- subset(clumpp.data, max) # subset and get best models
  setwd(wd)
  p <- ggplot(data = clumpp.data, aes(y=H, x=as.factor(model), fill=as.factor(perc.pca))) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values = colors, name="% Retained\nVariance") +
    ylab(label = "H'") +
    xlab(label = "Delimitation Model") +
    ylim(c(0,1)) +
    geom_text(aes(label = round(H, 2)),
              position = position_dodge(0.9),
              color="white",vjust = 1, hjust = 0.5,
              size=3) +
    theme_minimal()

  if(plot.type =="svg"){
    svglite(file = paste(plot.name, ".svg", sep=""), width = plot.width, height = plot.height)
    print(p)
    dev.off()
  } else if(plot.type =="pdf"){
    pdf(file = paste(plot.name, ".pdf", sep=""), width = plot.width, height = plot.height)
    print(p)
    dev.off()
  } else {
    png(file = paste(plot.name, ".png", sep=""), units="in", res=300, width = plot.width, height = plot.height)
    print(p)
    dev.off()
  }

  return(p)

}

