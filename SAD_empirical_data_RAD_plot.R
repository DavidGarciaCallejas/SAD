##################
# rank-abundance plots
##################

my.palette <- c("darkgreen","#009E73","#E69F00","#D55E00")

##################
##################

# abundances data
abundances <- readr::read_delim(file = "./data/abundance datasets/all_abundances.csv",delim = ";",col_types = list(col_character(),col_character(),col_double()))
abundances <- subset(abundances, !is.na(abundance))

# species data
species.data <- readr::read_delim(file = "./data/abundance datasets/all_species.csv",
                                  delim = ";",
                                  col_types = list(col_character(),
                                                   col_character(),
                                                   col_character(),
                                                   col_character(),
                                                   col_character()))

abundances$trophic.guild <- species.data$trophic.guild[match(abundances$species.ID,species.data$species.ID)]

# metrics data
all.metrics <- readr::read_delim(file = "./results/empirical_metrics.csv",delim = ";",col_types = list(col_character(),
                                                                                                           col_character(),
                                                                                                           col_double(),col_double(),col_double(),col_double(),
                                                                                                           col_double(),col_double(),col_double(),col_double(),
                                                                                                           col_double(),col_double(),col_double(),col_double()))

all.metrics$trophic.guild <- factor(all.metrics$trophic.guild, levels = c("plants","herbivores","omnivores","carnivores"))

site.parameters <- readr::read_delim(file = "./data/abundance datasets/all_sites.csv",delim = ";")

#- guilds with more than X species, otherwise it does not make sense to compute SAD
minimum.richness <- 3
all.metrics <- subset(all.metrics, richness >= minimum.richness)

# some cleaning
abundances <- subset(abundances,!is.na(trophic.guild))

abundances <- subset(abundances, abundance > 0) 
abundances$trophic.guild <- factor(abundances$trophic.guild, levels = c("plants","herbivores","omnivores","carnivores"))

##################
# rank-abundance curves
rank.abundances <- abundances[,c("site.ID","trophic.guild","abundance")]
rank.abundances <- rank.abundances %>% group_by(site.ID,trophic.guild) %>% mutate(relative.abund = abundance/sum(abundance))
rank.abundances <- rank.abundances %>% group_by(site.ID,trophic.guild) %>% mutate(species.rank = rank(-relative.abund,ties.method = "first"))

# all sites
# rad.plot <- ggplot(rank.abundances,aes(x = species.rank,y = relative.abund, color = trophic.guild, group = site.ID)) +
#   geom_point() + 
#   geom_line() +
#   facet_grid(.~trophic.guild,scales = "free")+
#   scale_color_manual(values = my.palette)+
#   theme_Publication()+
#   theme(strip.background = element_blank(),strip.text.x = element_blank()) +
#   guides(color=FALSE)+
#   NULL
# rad.plot

# select some sites
# trophic.guilds.per.site <- rank.abundances %>% count(site.ID,trophic.guild)
# trophic.guilds.per.site <- trophic.guilds.per.site %>% group_by(site.ID) %>% summarise(num.guilds = n())
# diverse.sites <- trophic.guilds.per.site$site.ID[trophic.guilds.per.site$num.guilds > 2]
# 
# # are the guilds well balanced?
# diverse.sites.abund <- subset(abundances,site.ID %in% diverse.sites)
# diverse.sites.abund <- diverse.sites.abund %>% group_by(site.ID,trophic.guild) %>% summarise(num.sp = n())
# 
# diverse.sites.plot <- ggplot(diverse.sites.abund, aes(x = trophic.guild, y = num.sp)) +
#   geom_point() + 
#   facet_wrap(~site.ID) +
#   NULL
# diverse.sites.plot
# 
# # sites with all guilds with >1 species
# one.sp.sites <- unique(diverse.sites.abund$site.ID[diverse.sites.abund$num.sp == 1])
# diverse.sites.abund <- subset(diverse.sites.abund,!(site.ID %in% one.sp.sites))
# 
# diverse.sites.plot <- ggplot(diverse.sites.abund, aes(x = trophic.guild, y = num.sp)) +
#   geom_point() + 
#   facet_wrap(~site.ID) + ylim(c(3,20)) +
#   NULL
# diverse.sites.plot

selected.sites <- c("ALTERDOC",1200,1271,1282,1822)

selected.abundances <- subset(rank.abundances,site.ID %in% selected.sites)
selected.abundances$site.ID <- factor(selected.abundances$site.ID, levels = c("ALTERDOC",1200,1271,1282,1822))

selected.rad.plot <- ggplot(selected.abundances,aes(x = species.rank,y = relative.abund, group = trophic.guild)) +
  geom_line(aes(color = trophic.guild), size = 1.1) +#(aes(linetype = site.ID)) +
  geom_point(aes(fill = trophic.guild),shape = 21, size = 1.5) +#(aes(shape = site.ID)) + 
  facet_grid(.~site.ID,scales = "free")+
  scale_color_manual(values = my.palette)+#, labels = c("plants", "herbivores", "omnivores", "carnivores (inv)", "carnivores (vert)"))+
  scale_fill_manual(values = my.palette)+#, labels = c("plants", "foli/granivores", "frugi/nectarivores", "omnivores", "carnivores (inv)", "carnivores (vert)"))+
                      #c("primary\nproducers", "plant/seed\neaters", "frugivores", "omnivores", "carnivores\n(invertebrates)", "carnivores\n(vertebrates)"))+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
  xlab("species rank") + ylab("relative abundance") +
  #DGC::theme_Publication()+
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  #guides(color=FALSE)+#, fill = FALSE)+
  NULL

tiff("./results/images/empirical_RAD_sites.tiff", res=600, compression = "lzw", width = 5500, height = 2500, units = "px")
selected.rad.plot
dev.off()

###########################################
###########################################
# here I draw rad plots of highly even and uneven sites with at least four species

# all.metrics <- droplevels(subset(all.metrics, richness > 4 & hill.evenness < 0.99))
# 
# evenness.plants <- all.metrics$hill.evenness[all.metrics$trophic.guild == "plants"]
# evenness.plants <- sort(evenness.plants,decreasing = T)
# 
# plants.highest <- all.metrics$site.ID[all.metrics$trophic.guild == "plants" & all.metrics$hill.evenness %in% evenness.plants[1:2]]
# plants.lowest <- all.metrics$site.ID[all.metrics$trophic.guild == "plants" & all.metrics$hill.evenness %in% evenness.plants[(length(evenness.plants)-1):length(evenness.plants)]]
# 
# evenness.herbivores <- all.metrics$hill.evenness[all.metrics$trophic.guild == "herbivores"]
# evenness.herbivores <- sort(evenness.herbivores,decreasing = T)
# 
# herbivores.highest <- all.metrics$site.ID[all.metrics$trophic.guild == "herbivores" & all.metrics$hill.evenness %in% evenness.herbivores[1:2]]
# herbivores.lowest <- all.metrics$site.ID[all.metrics$trophic.guild == "herbivores" & all.metrics$hill.evenness %in% evenness.herbivores[(length(evenness.herbivores)-1):length(evenness.herbivores)]]
# 
# evenness.omnivores <- all.metrics$hill.evenness[all.metrics$trophic.guild == "omnivores"]
# evenness.omnivores <- sort(evenness.omnivores,decreasing = T)
# 
# omnivores.highest <- all.metrics$site.ID[all.metrics$trophic.guild == "omnivores" & all.metrics$hill.evenness %in% evenness.omnivores[1:2]]
# omnivores.lowest <- all.metrics$site.ID[all.metrics$trophic.guild == "omnivores" & all.metrics$hill.evenness %in% evenness.omnivores[(length(evenness.omnivores)-1):length(evenness.omnivores)]]
# 
# evenness.carnivores <- all.metrics$hill.evenness[all.metrics$trophic.guild == "carnivores"]
# evenness.carnivores <- sort(evenness.carnivores,decreasing = T)
# 
# carnivores.highest <- all.metrics$site.ID[all.metrics$trophic.guild == "carnivores" & all.metrics$hill.evenness %in% evenness.carnivores[1:2]]
# carnivores.lowest <- all.metrics$site.ID[all.metrics$trophic.guild == "carnivores" & all.metrics$hill.evenness %in% evenness.carnivores[(length(evenness.carnivores)-1):length(evenness.carnivores)]]
# 
# ###############
# plants.rad <- abundances[abundances$site.ID %in% c(plants.highest,plants.lowest),]
# herb.rad <- abundances[abundances$site.ID %in% c(herbivores.highest,herbivores.lowest) & abundances$trophic.guild == "herbivores",]
# omni.rad <- abundances[abundances$site.ID %in% c(omnivores.highest,omnivores.lowest) & abundances$trophic.guild == "omnivores",]
# carn.rad <- abundances[abundances$site.ID %in% c(carnivores.highest,carnivores.lowest) & abundances$trophic.guild == "carnivores",]
# evenness.rad <- bind_rows(plants.rad,herb.rad,omni.rad,carn.rad)
# 
# evenness.rad$evenness <- "low"
# evenness.rad$evenness[evenness.rad$site.ID %in% plants.highest & evenness.rad$trophic.guild == "plants"] <- "high"
# evenness.rad$evenness[evenness.rad$site.ID %in% herbivores.highest & evenness.rad$trophic.guild == "herbivores"] <- "high"
# evenness.rad$evenness[evenness.rad$site.ID %in% omnivores.highest & evenness.rad$trophic.guild == "omnivores"] <- "high"
# evenness.rad$evenness[evenness.rad$site.ID %in% carnivores.highest & evenness.rad$trophic.guild == "carnivores"] <- "high"
# 
# evenness.rad <- evenness.rad[,c("site.ID","trophic.guild","abundance","evenness")]
# evenness.rad <- evenness.rad %>% group_by(site.ID,trophic.guild) %>% mutate(relative.abund = abundance/sum(abundance))
# evenness.rad <- evenness.rad %>% group_by(site.ID,trophic.guild) %>% mutate(species.rank = rank(-relative.abund,ties.method = "first"))
# ###############
# 
# evenness.rad.plot <- ggplot(evenness.rad,aes(x = species.rank,y = relative.abund, group = site.ID)) +
#   geom_line(aes(color = trophic.guild), size = 1.1) +#(aes(linetype = site.ID)) +
#   geom_point(aes(fill = trophic.guild),shape = 21, size = 1.5) +#(aes(shape = site.ID)) + 
#   facet_wrap(trophic.guild~evenness,scales = "free")+
#   scale_color_manual(values = my.palette)+#, labels = c("plants", "herbivores", "omnivores", "carnivores (inv)", "carnivores (vert)"))+
#   scale_fill_manual(values = my.palette)+#, labels = c("plants", "foli/granivores", "frugi/nectarivores", "omnivores", "carnivores (inv)", "carnivores (vert)"))+
#   #c("primary\nproducers", "plant/seed\neaters", "frugivores", "omnivores", "carnivores\n(invertebrates)", "carnivores\n(vertebrates)"))+
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
#   xlab("species rank") + ylab("relative abundance") +
#   #DGC::theme_Publication()+
#   theme(strip.background = element_blank())+#,strip.text.x = element_blank()) +
#   guides(color=FALSE, fill = FALSE)+
#   NULL
# # evenness.rad.plot
