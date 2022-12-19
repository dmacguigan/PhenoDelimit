# plot permuted H' vs observed H'
H_permutation_plot <- function(wd, clumpp_results, clumpp_perm_df, models, best.perc.var, plot.type, plot.prefix, plot.width, plot.height, sig.threshold){
  setwd(wd)
  summary <- data.frame(model=models,
                        H_diff=models,
                        p_val=models)
  counter=1
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
    
    # calculate p value
    p_val <- (sum(perm_H$H > obs_H$H) + 1)/nrow(perm_H) # p-value is the proportion of samples (including our observed data) that have a test statistic larger than that of our observed data 
    message(p_val)
    
    # calculate diff between obs H' and mean permuted H'
    H_diff <- obs_H$H - mean(perm_H$H)
    message(obs_H$H)
    message(mean(perm_H$H))
    message(H_diff)
    
    # update the summary data frame
    summary[counter,2] = H_diff
    summary[counter,3] = p_val
    
    counter=counter+1
  }
  
  # plot CLUMPP H' obs - permutation mean for each delimitation model
  p <- ggplot(data = summary, aes(y=H_diff, x=as.factor(model))) +
    geom_bar(stat="identity", position=position_dodge(), fill="black") +
    ylab(label = "Observed H' - Permutation Mean H'") +
    xlab(label = "Delimitation Model") +
    geom_text(aes(label = round(H_diff, 2)),
              position = position_dodge(0.9),
              color="white",vjust = 1.1, hjust = 0.5,
              size=3) +
    geom_text(aes(label = round(p_val, 3)),
              color="black",vjust = -0.5, hjust = 0.5,
              size=3) +
    theme_minimal()
  
  plot.name <- paste0(plot.prefix, "_obs-minus-perm-mean")
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
  
}
