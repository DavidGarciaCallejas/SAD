###################
# plots from the theoretical model
###################

my.palette <- c("darkgreen","#009E73","#E69F00","#D55E00")

# metrics.results <- readr::read_delim(file = "./results/sugihara_model_metrics_DD.csv",delim = ";")
# metrics.results <- readr::read_delim(file = "./results/sugihara_model_metrics_DP.csv",delim = ";")
# metrics.results <- readr::read_delim(file = "./results/sugihara_model_metrics_RF.csv",delim = ";")
metrics.results <- readr::read_delim(file = "./results/sugihara_model_metrics_all_apportionments.csv",delim = ";")


richness.levels <- unique(metrics.results$richness.level)
connectance.levels <- unique(metrics.results$connectance.level)
apportionment.levels <- unique(metrics.results$apportionment.level)
replicates <- max(metrics.results$replicate)
trophic.guilds <- unique(metrics.results$trophic.guild)

metrics.results$richness.level <- factor(metrics.results$richness.level)
metrics.results$connectance.level <- factor(metrics.results$connectance.level)
metrics.results$apportionment.level <- factor(metrics.results$apportionment.level)
metrics.results$trophic.guild <- factor(metrics.results$trophic.guild, levels = c("primary Producers","herbivores","omnivores","carnivores"))

#####################
### APPLY any final subset and/or renaming here ###
# you will have to modify the plot calling, in particular the facetting, according to the subsetting here

# for these analyses, I focus on random fraction apportionment
plot.data <- metrics.results
# plot.data <- subset(metrics.results,apportionment.level == "random.fraction" & richness.level == 100 & connectance.level == 0.2)

plot.data$trophic.guild <- plyr::revalue(plot.data$trophic.guild,c("primary Producers" = "primary\n producers"))

#####################

metric.centroids <- plot.data %>% 
  # dplyr::filter(apportionment.level == apportionment.levels[i.apport]) %>% 
  group_by(richness.level,
           apportionment.level,
           connectance.level,
           trophic.guild) %>% 
  summarize(centroid.skewness = mean(skewness),centroid.evenness = mean(hill.evenness))

#####################
# plot.list <- list()
# 
# for(i.apport in 1:length(apportionment.levels)){

sugihara.metrics.plot <- ggplot(plot.data,aes(x = skewness,y = hill.evenness)) +
  
  geom_vline(xintercept = 0,color="grey80") +
  geom_hline(yintercept = 1,color="grey80") +
  
  geom_point(aes(color = trophic.guild), size = 1.2, alpha = 0.7) +
  geom_point(data = metric.centroids,aes(x = centroid.skewness, y = centroid.evenness, fill = trophic.guild), size = 3, shape = 21) +
  
  # facet_grid(richness.level~trophic.guild+connectance.level) +
  facet_grid(connectance.level~richness.level+trophic.guild) +
  # facet_grid(.~trophic.guild) +
  
  scale_fill_manual(values = my.palette)+
  scale_color_manual(values = my.palette)+
  theme(panel.grid.minor=element_blank()) +
  # DGC::theme_Publication() +
  theme(strip.background = element_blank()) +
  # theme(panel.border = element_rect(colour="black")) +
  guides(fill=F,color = F) +
  
  xlim(-1,1) + #ylim(.45,1) + 
  xlab("skewness") + ylab("evenness") +
  # ggtitle("",subtitle = apportionment.levels[i.apport]) +
  NULL

# plot.list[[i.apport]] <- sugihara.metrics.plot

# }

# tiff("./results/images/sugihara_metrics_DD.tiff", res=600, compression = "lzw", width = 4500, height = 2000, units = "px")
# tiff("./results/images/sugihara_metrics_DP.tiff", res=600, compression = "lzw", width = 4500, height = 2000, units = "px")
# tiff("./results/images/sugihara_metrics_RF.tiff", res=600, compression = "lzw", width = 4500, height = 2000, units = "px")

sugihara.metrics.plot
# dev.off()
