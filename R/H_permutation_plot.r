# plot permuted H' vs observed H'
H_permutation_plot <- function(clumpp_results, clumpp_perm_df, models, best.perc.var, plot.type, plot.prefix, plot.width, plot.height){
  for(i in models){
    perm_H <- subset(clumpp_perm_df, model==i & perc.pca==best.perc.var)
    obs_H <- subset(clumpp_results, model==i & perc.pca==best.perc.var)
    p <- ggplot(data=perm_H, aes(x=H)) +
      geom_histogram(color="black", fill="gray") +
      xlab("H'") +
      ggtitle(paste0("model ", i, ", ", best.perc.var, "% retained variance")) +
      geom_vline(data=obs_H, aes(xintercept=H), color="red", linetype="dashed", size=1.5)+
      theme_minimal()
    
    plot.name <- paste0(plot.prefix, "_perm-vs-obs_m", i, "")
    
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
      
    #um(dist > diff(by(y, tr, mean)))/2000  #
  }
}