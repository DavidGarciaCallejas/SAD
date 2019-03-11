###################
# plots from the theoretical model
###################

# metrics.results <- readr::read_delim(file = "./results/sugihara_model_metrics_DD.csv",delim = ";")
# metrics.results <- readr::read_delim(file = "./results/sugihara_model_metrics_DP.csv",delim = ";")
# metrics.results <- readr::read_delim(file = "./results/sugihara_model_metrics_RF.csv",delim = ";")
metrics.results <- readr::read_delim(file = "./results/sugihara_model_metrics_complete.csv",delim = ";")

resource.dist.levels <- unique(metrics.results$resource.distribution.level)
richness.levels <- unique(metrics.results$richness.level)
connectance.levels <- unique(metrics.results$connectance.level)
apportionment.levels <- unique(metrics.results$apportionment.level)
replicates <- max(metrics.results$replicate)
trophic.guilds <- unique(metrics.results$trophic.guild)

metrics.results$resource.distribution.level <- factor(metrics.results$resource.distribution.level,levels = c("uniform","skewed"))
metrics.results$richness.level <- factor(metrics.results$richness.level)
metrics.results$connectance.level <- factor(metrics.results$connectance.level, levels = c("0.3","0.2","0.1"))
metrics.results$apportionment.level <- recode(metrics.results$apportionment.level,"random.assortment" = "random assortment",
                                              "random.fraction" = "random fraction",
                                              "dominance.decay" = "dominance decay",
                                              "dominance.preemption" = "dominance preemption")
# reorder levels
metrics.results$apportionment.level <- factor(metrics.results$apportionment.level,levels = c("random assortment","dominance decay","random fraction","dominance preemption"))
# metrics.results$trophic.guild <- factor(metrics.results$trophic.guild, levels = c("primary Producers","herbivores","omnivores","carnivores"))
metrics.results$trophic.guild <- factor(metrics.results$trophic.guild,levels = c("resources","basal","intermediate","top"))

#####################
### APPLY any final subset and/or renaming here ###
# you will have to modify the plot calling, in particular the facetting, according to the subsetting here

# for these analyses, I focus on random fraction apportionment
# plot.data <- metrics.results
plot.data <- subset(metrics.results,resource.distribution.level == "skewed" & trophic.guild != "resources")

#####################

# metric.centroids <- plot.data %>% 
#   # dplyr::filter(apportionment.level == apportionment.levels[i.apport]) %>% 
#   group_by(richness.level,
#            apportionment.level,
#            connectance.level,
#            trophic.guild) %>% 
#   summarize(centroid.skewness = mean(skewness),centroid.evenness = mean(hill.evenness))

#####################
# plot.list <- list()
# 
# for(i.apport in 1:length(apportionment.levels)){

my.palette <- c("gray60","#009E73","#E69F00","#D55E00")
# my.palette <- c("gray70","gray50","gray30","gray10")
# library(ggbeeswarm)

mean.data <- plot.data %>% group_by(trophic.guild,apportionment.level,connectance.level) %>% summarise(mean.ev = mean(hill.evenness),sd.ev = sd(hill.evenness))
pos = position_dodge(.2)

evenness.mean.plot <- ggplot(mean.data,aes(x = trophic.guild, y = mean.ev,group = apportionment.level,colour = apportionment.level, fill = apportionment.level))+
  geom_line(position = pos, size = 1.1)+
  geom_errorbar(aes(ymin = mean.ev - sd.ev, ymax = mean.ev+sd.ev),width = .3,position = pos) +
  geom_point(position = pos, aes(shape = apportionment.level), color = "black",size = 1.7) + 
  
  facet_grid(connectance.level~.)+
  scale_shape_manual(values = c(21,22,23,24))+
  scale_fill_manual(values = my.palette)+
  scale_color_manual(values = my.palette)+
  DGC::theme_Publication()+
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(legend.position="none")+
  xlab("trophic guild") + ylab("evenness")+
  NULL
# evenness.mean.plot

# sugihara.evenness.plot <- ggplot(plot.data, aes(x = trophic.guild, y = hill.evenness, group = apportionment.level,colour = apportionment.level, fill = apportionment.level))+
#   # stat_summary(fun.data = "mean_se", geom="errorbar", width = .05) + 
#   stat_summary(fun.y = mean,
#                fun.ymin = function(x) mean(x) - sd(x),
#                fun.ymax = function(x) mean(x) + sd(x),
#                geom = "pointrange",
#                position = position_dodge(.1),size = .5,aes(shape = apportionment.level),color = "black") +
#   stat_summary(fun.y = mean,
#                geom = "line", position = position_dodge(.1), size = 1.2) +
#   # geom_boxplot(aes(fill = apportionment.level),outlier.size = .5) +
#   # geom_violin(aes(fill = apportionment.level))+
#   # geom_quasirandom(aes(fill = apportionment.level)) +
#   facet_grid(connectance.level~.)+
#   scale_shape_manual(values = c(21,22,23,24))+
#   scale_fill_manual(values = my.palette)+
#   scale_color_manual(values = my.palette)+
#   DGC::theme_Publication()+
#   theme(strip.background = element_blank(),strip.text.x = element_blank()) +
#   xlab("trophic guild") + ylab("evenness")+
#   NULL
# 
# sugihara.evenness.plot

# sugihara.evenness.plot <- ggplot(plot.data, aes(x = trophic.guild, y = hill.evenness))+
#   geom_boxplot(aes(fill = apportionment.level),outlier.size = .5) +
#   # geom_violin(aes(fill = apportionment.level))+
#   # geom_quasirandom(aes(fill = apportionment.level)) +
#   facet_grid(connectance.level~.)+
#   scale_fill_manual(values = my.palette)+
#   DGC::theme_Publication()+
#   theme(strip.background = element_blank(),strip.text.x = element_blank()) +
#   theme(legend.position="none")+
#   xlab("trophic guild") + ylab("evenness")+
#   NULL


tiff("./results/images/summary_evenness.tiff", res=600, compression = "lzw", width = 3500, height = 3000, units = "px")
evenness.mean.plot
dev.off()

