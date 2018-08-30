##################
# plot metrics from empirical datasets
##################

my.palette <- c("darkgreen","#009E73","#E69F00","#D55E00")

##################
##################

# metrics data
all.metrics <- readr::read_delim(file = "./results/empirical_metrics.csv",delim = ";",col_types = list(col_character(),
                                                                                                      col_character(),
                                                                                                      col_double(),col_double(),col_double(),col_double(),
                                                                                                      col_double(),col_double(),col_double(),col_double(),
                                                                                                      col_double(),col_double(),col_double(),col_double()))

all.metrics$trophic.guild <- factor(all.metrics$trophic.guild, levels = c("plants","herbivores","omnivores","carnivores"))
# 
# all.metrics$trophic.guild <- plyr::revalue(all.metrics$trophic.guild, c("PrimaryProducer"="plants", 
#                                                                                     "PlantSeed"="foli/granivores",
#                                                                                     "FruiNect" = "frugi/nectarivores",
#                                                                                     "Omnivore" = "omnivores",
#                                                                                     "Invertebrate" = "carnivores (inv)",
#                                                                                     "VertFishScav" = "carnivores (vert)"))

site.parameters <- readr::read_delim(file = "./data/abundance datasets/all_sites.csv",delim = ";")

#- guilds with more than X species, otherwise it does not make sense to compute SAD
minimum.richness <- 3
all.metrics <- subset(all.metrics, richness >= minimum.richness)
############
# natural or log-transformed abundances?
# empirical.metrics <- empirical.metrics[,c("site.ID","trophic.guild","hill.evenness.log","skewness.log")]
# names(empirical.metrics)[c(3,4)] <- c("hill.evenness","skewness")
############
empirical.metrics <- all.metrics[,c("site.ID","trophic.guild","hill.evenness","skewness")]
############


empirical.metrics <- gather(empirical.metrics,key = "metric",value = "value",-site.ID,-trophic.guild)

# plot separately for better tuning of the plots
evenness.data <- subset(empirical.metrics,metric == "hill.evenness")
skewness.data <- subset(empirical.metrics,metric == "skewness")

##############
evenness.plot <- ggplot(evenness.data) +
  geom_density_ridges(aes(x = value, y = trophic.guild, fill = trophic.guild, point_fill = trophic.guild),
                      point_shape = 21,
                      point_size = 0.7,
                      alpha = .7, jittered_points = TRUE,
                      # rel_min_height = 0.001,
                      scale = 0.9) +
  
  scale_fill_manual(values = my.palette)+
  scale_discrete_manual(aesthetics = "point_fill", values = my.palette) +
  # DGC::theme_Publication()+
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  guides(fill=FALSE,point_fill = FALSE)+
  xlab("evenness")+ ylab("") +
  scale_y_discrete(breaks=NULL,expand = c(0.01, 0)) +
  scale_x_continuous(limits=c(0,1),expand = c(0.01, 0))+ # when using ggridges
  NULL

# tiff("/.results/images/empirical_evenness.tiff", res=600, compression = "lzw", width = 5500, height = 2500, units = "px")
#  evenness.plot
# dev.off()

################

skewness.plot <- ggplot(skewness.data) +
   geom_density_ridges(aes(x = value, y = trophic.guild, fill = trophic.guild, point_fill = trophic.guild),
                       point_shape = 21,
                       point_size = 0.7,
                       alpha = .7, jittered_points = TRUE,
                       # rel_min_height = 0.001,
                       scale = 0.9) +
  scale_fill_manual(values = my.palette)+
  scale_discrete_manual(aesthetics = "point_fill", values = my.palette) +
  # DGC::theme_Publication()+
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  guides(fill=FALSE,point_fill = FALSE)+
  xlab("skewness")+ ylab("") + 
  scale_y_discrete(breaks=NULL,expand = c(0.01, 0)) +
  scale_x_continuous(limits=c(-1,1),expand = c(0.01, 0))+ # when using ggridges
  NULL

# tiff("./results/images/empirical_skewness.tiff", res=600, compression = "lzw", width = 5500, height = 2500, units = "px")
# skewness.plot
# dev.off()

tiff("./results/images/empirical_metrics_ridge.tiff", res=600, compression = "lzw", width = 4500, height = 2500, units = "px")
evenness.plot + skewness.plot
dev.off()

#############
# how does it look when plotted as the cohen results?

data.centroids <- empirical.metrics %>% group_by(trophic.guild,metric) %>% summarize(centroid = mean(value))
data.centroids <- data.centroids %>% spread(key = metric,value = centroid)
empirical.metrics <- empirical.metrics %>% spread(metric,value)

empirical.metrics.plot <- ggplot(empirical.metrics,aes(x = skewness,y = hill.evenness, fill = trophic.guild))
empirical.metrics.plot <- empirical.metrics.plot +
  
  geom_vline(xintercept = 0,color="grey80") +
  geom_hline(yintercept = 1,color="grey80") +
  
  geom_point(size = 1.2, alpha = .5,shape = 21) +
  geom_point(data = data.centroids,aes(x = skewness, y = hill.evenness), fill = "grey80", size = 3, shape = 21) +
  
  scale_fill_manual(values = my.palette) +
  
  facet_grid(.~trophic.guild) +
  # theme(panel.grid.minor=element_blank()) +
  # DGC::theme_Publication() + 
  theme(strip.background = element_blank()) +
  # theme(panel.border = element_rect(colour="black")) +
  guides(fill=F) +
  
  xlim(-1,1) + #ylim(.45,1) + 
  xlab("skewness") + ylab("evenness") +
  NULL
# empirical.metrics.plot

tiff("./results/images/empirical_metrics_complete.tiff", res=600, compression = "lzw", width = 4500, height = 3500, units = "px")
evenness.plot + skewness.plot - empirical.metrics.plot + plot_layout(ncol = 1)
dev.off()


