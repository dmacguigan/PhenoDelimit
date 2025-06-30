#' Plot CLUMPP H' values
#' 
#' Plot CLUMPP H' values for each delimitation model and percetage of retained variance. Returns plot and saves plot as file (svg or pdf).
#' 
#' @param wd directory to create plots
#' @param clumpp.data data frame from read_clumpp_results function
#' @param colors vector of colors for percentage of retained variance
#' @param plot.name name for plot, character
#' @param plot.type save plot as "pdf", "svg", or "png"
#' @param plot.width width of plot in inches
#' @param plot.height height of plot in inches
#' 
#' @export
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

