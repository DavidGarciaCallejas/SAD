###################
# plots from the theoretical model
###################

# metrics.results <- readr::read_delim(file = "./results/cohen_simulation_metrics_CONSTANT_DEGREE.csv",delim = ";")
metrics.results <- readr::read_delim(file = "./results/cohen_simulation_metrics_VARIABLE_DEGREE.csv",delim = ";")

vulnerability.levels <- unique(metrics.results$vulnerability.level)
bodysize.tl.levels <- unique(metrics.results$body.size.tl)
replicates <- max(metrics.results$replicate)
trophic.levels <- max(metrics.results$trophic.level)

metrics.results$vulnerability.level <- factor(metrics.results$vulnerability.level,levels = c("min","low","intermediate","high","max"))
metrics.results$body.size.tl <- factor(metrics.results$body.size.tl)
# metrics.results$body.size.tl <- factor(metrics.results$body.size.tl,levels = c("decreasing","none","increasing"))

#####################
### APPLY any final subset and renaming here ###
# plot.data <- metrics.results
# plot.data <- subset(plot.data,vulnerability.level %in% c("min","low","intermediate"))

plot.data <- droplevels(subset(metrics.results,vulnerability.level %in% c("low","intermediate","high")))

cohen.centroids <- plot.data %>% group_by(vulnerability.level,body.size.tl,trophic.level) %>% summarize(centroid.skewness = mean(skewness),centroid.evenness = mean(hill.evenness))
centroids.plot <- cohen.centroids

#####################
#####################

cohen.metrics.plot <- ggplot(plot.data,aes(x = skewness,y = hill.evenness)) +

  geom_vline(xintercept = 0,color="grey80") +
  geom_hline(yintercept = 1,color="grey80") +
  
  geom_point(size = 0.6, alpha = 0.4) +
  geom_point(data = centroids.plot,aes(x = centroid.skewness, y = centroid.evenness), fill = "#E69F00", size = 3, shape = 21) +
  
  facet_grid(vulnerability.level~trophic.level) +
  # theme(panel.grid.minor=element_blank()) +
  # DGC::theme_Publication() + 
  theme(strip.background = element_blank()) +
  # theme(panel.border = element_rect(colour="black")) +
  # guides(fill=F) +
  
  xlim(-1,1) + #ylim(.45,1) + 
  xlab("skewness") + ylab("evenness") +
  NULL

tiff("./results/images/cohen_metrics_VARIABLE_DEGREE.tiff", res=600, compression = "lzw", width = 4500, height = 4000, units = "px")
cohen.metrics.plot
dev.off()
