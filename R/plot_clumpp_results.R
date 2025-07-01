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
  p <- ggplot2::ggplot(data = clumpp.data, ggplot2::aes(y=H, x=as.factor(model), fill=as.factor(perc.pca))) +
    ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge()) +
    ggplot2::scale_fill_manual(values = colors, name="% Retained\nVariance") +
    ggplot2::ylab(label = "H'") +
    ggplot2::xlab(label = "Delimitation Model") +
    ggplot2::ylim(c(0,1)) +
    ggplot2::geom_text(ggplot2::aes(label = round(H, 2)),
              position = ggplot2::position_dodge(0.9),
              color="white",vjust = 1, hjust = 0.5,
              size=3) +
    ggplot2::theme_minimal()

  if(plot.type =="svg"){
    svglite::svglite(file = paste(plot.name, ".svg", sep=""), width = plot.width, height = plot.height)
    print(p)
    dev.off()
  } else if(plot.type =="pdf"){
    grDevices::pdf(file = paste(plot.name, ".pdf", sep=""), width = plot.width, height = plot.height)
    print(p)
    dev.off()
  } else {
    grDevices::png(file = paste(plot.name, ".png", sep=""), units="in", res=300, width = plot.width, height = plot.height)
    print(p)
    dev.off()
  }

  return(p)

}

