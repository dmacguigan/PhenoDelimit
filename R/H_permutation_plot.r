#' Plot permuted CLUMPP results
#' 
#' Plot permuted CLUMPP H' results. Returns plot of obsserved H' minus mean permuted H'. Number inside bar is the H' difference, number on top of bar is p-value from permutation test. 
#` Also returns plots for each model of permuted H' distribution and observed H' (red dashed line). Saves plot as file (svg or pdf).
#' 
#' @param wd directory to create plots, should contain "CLUMPP_permuted" subdirectory
#' @param clumpp.data data frame from read_clumpp_results function
#' @param clumpp.data.permuted data frame from read_clumpp_results_permuted function
#' @param model.numbers vector containing delimitaiton model numbers
#' @param best.perc.var which percent retained variance to plot? numeric
#' @param plot.prefix prefix for plot name, character
#' @param plot.type save plot as "pdf", "svg", or "png"
#' @param plot.width width of plot in inches
#' @param plot.height height of plot in inches
#' 
#' @export

H_permutation_plot <- function(wd, clumpp.data, clumpp.data.permuted, model.numbers, best.perc.var, plot.type, plot.prefix, plot.width, plot.height){
  setwd(wd)
  summary <- data.frame(model=model.numbers,
                        H_diff=model.numbers,
                        p_val=model.numbers)
  counter=1
  for(i in model.numbers){
    perm_H <- subset(clumpp.data.permuted, model==i & perc.pca==best.perc.var)
    obs_H <- subset(clumpp.data, model==i & perc.pca==best.perc.var)
    p <- ggplot2::ggplot(data=perm_H, ggplot2::aes(x=H)) +
      ggplot2::geom_histogram(color="black", fill="gray") +
      ggplot2::xlab("H'") +
      ggplot2::ggtitle(paste0("model ", i, ", ", best.perc.var, "% retained variance")) +
      ggplot2::geom_vline(data=obs_H, ggplot2::aes(xintercept=H), color="red", linetype="dashed", linewidth=1.5)+
      ggplot2::theme_minimal()
    
    plot.name <- paste0(plot.prefix, "_perm-vs-obs_m", i, "")
    
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
    
    # calculate p value
    p_val <- (sum(perm_H$H > obs_H$H) + 1)/nrow(perm_H) # p-value is the proportion of samples (including our observed data) that have a test statistic larger than that of our observed data 
    
    # calculate diff between obs H' and mean permuted H'
    H_diff <- obs_H$H - mean(perm_H$H)
    
    # update the summary data frame
    summary[counter,2] = H_diff
    summary[counter,3] = p_val
    
    counter=counter+1
  }
  
  # plot CLUMPP H' obs - permutation mean for each delimitation model
  p <- ggplot2::ggplot(data = summary, ggplot2::aes(y=H_diff, x=as.factor(model))) +
    ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge(), fill="black") +
    ggplot2::ylab(label = "Observed H' - Permutation Mean H'") +
    ggplot2::xlab(label = "Delimitation Model") +
    ggplot2::geom_text(ggplot2::aes(label = round(H_diff, 2)),
              position = ggplot2::position_dodge(0.9),
              color="white",vjust = 1.1, hjust = 0.5,
              size=3) +
    ggplot2::geom_text(ggplot2::aes(label = paste0("p=",round(p_val, 3))),
              color="black",vjust = -0.5, hjust = 0.5,
              size=3) +
    ggplot2::theme_minimal()
  
  plot.name <- paste0(plot.prefix, "_obs-minus-perm-mean")
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
  
}
