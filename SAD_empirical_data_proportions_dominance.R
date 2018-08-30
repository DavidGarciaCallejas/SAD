##################
# rank-abundance plots
##################

# my.palette <- c("darkgreen","#009E73","#E69F00","#D55E00") # no plants in these analyses
my.palette <- c("#009E73","#E69F00","#D55E00")

##################

p_dodge <- position_dodge(2)

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

############
# calculate indices for every trophic guild of each site
abundances$trophic.guild <- species.data$trophic.guild[match(abundances$species.ID,species.data$species.ID)]

# some cleaning
abundances <- subset(abundances,!is.na(trophic.guild))
abundances <- subset(abundances, abundance > 0 & trophic.guild != "plants")
abundances$trophic.guild <- factor(abundances$trophic.guild, levels = c("herbivores","omnivores","carnivores"))
###############
# individuals and richness proportions
site.proportions <- abundances %>% group_by(site.ID,trophic.guild) %>% summarise(guild.abundance = sum(abundance),guild.richness = n())
site.proportions <- site.proportions %>% group_by(site.ID) %>% mutate(net.abundance = sum(guild.abundance),net.richness = sum(guild.richness))
site.proportions.stat <- site.proportions
site.proportions <- site.proportions %>% group_by(site.ID,trophic.guild) %>% summarise(abundance.proportion = guild.abundance/net.abundance, richness.proportion = guild.richness/net.richness)
mean.proportions <- site.proportions %>% group_by(trophic.guild) %>% summarise(mean.sp = mean(richness.proportion), mean.ind = mean(abundance.proportion))

# dominance
dominance.data <- abundances %>% group_by(site.ID,trophic.guild) %>% mutate(guild.abundance = sum(abundance))
dominance.data$sp.dominance <- dominance.data$abundance/dominance.data$guild.abundance
dominance.data <- dominance.data %>% group_by(site.ID,trophic.guild) %>% mutate(dominance = max(sp.dominance),guild.richness = n())
dominance.data <- dominance.data %>% group_by(trophic.guild,guild.richness) %>% mutate(mean.dominance = mean(dominance))

################
# plots
proportions.plot <- ggplot(site.proportions,aes(x = richness.proportion,y=abundance.proportion)) +
  geom_point(aes(fill = trophic.guild),shape = 21, alpha = 0.5, size = 1.6) +
  geom_abline(slope = 1,intercept = 0, color = "grey20") +
  geom_point(data = mean.proportions, aes(x = mean.sp, y = mean.ind), color = "black", fill = "grey80", shape = 21, size = 3)+
  facet_grid(.~trophic.guild)+
  #DGC::theme_Publication() +
  xlab("species proportion") + ylab("individuals proportion") +
  scale_fill_manual(values = my.palette)+
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  guides(fill = FALSE) +
  NULL


tiff("./results/images/ind_sp_proportions.tiff", res=600, compression = "lzw", width = 5500, height = 2500, units = "px")
proportions.plot
dev.off()

###############
# dominance
dominance.plot <- ggplot(dominance.data,aes(x = guild.richness,y = dominance,fill = trophic.guild)) +
  # geom_point(shape = 21, size = 0.7, alpha = .5) +
  geom_jitter(aes(fill = trophic.guild),position = position_jitter(width = 0.3, height = 0),shape = 21,size = 1.5,alpha =0.5,color = "black")+
  xlim(c(1,20))+ ylim(c(0,1))+

  facet_grid(.~trophic.guild)+
  stat_function(fun = function(x) 1/x, color = "grey20") +
  geom_hline(yintercept = 1, color = "grey20") +
  geom_point(aes(x = guild.richness, y = mean.dominance), color = "black", fill = "grey80", shape = 21, size = 3)+

  scale_fill_manual(values = my.palette)+
  # scale_color_manual(values = my.palette)+
  xlab("guild richness") + ylab("dominance") +
  #DGC::theme_Publication() +
  theme(strip.background = element_blank())+#,strip.text.x = element_blank()) +
  guides(fill = FALSE) +
  NULL

tiff("./results/images/dominance_richness.tiff", res=600, compression = "lzw", width = 5500, height = 2500, units = "px")
dominance.plot
dev.off()

# as an alternative, save the two plots together with patchwork, as in the paper:  dominance.plot + proportions.plot + plot_layout(ncol = 1)

#########################
# statistical analyses
# 1 - ind/sp proportions

site.proportions.stat <- site.proportions.stat %>% group_by(site.ID) %>% mutate(net.abundance = sum(guild.abundance),net.richness = sum(guild.richness))
site.proportions.stat <- site.proportions.stat %>% group_by(site.ID, trophic.guild) %>% summarise(species.proportion = guild.richness/net.richness, 
                                                                                                  individuals.proportion = guild.abundance/net.abundance)
proportions.test <- site.proportions.stat %>% group_by(trophic.guild) %>% do(w = wilcox.test(.$individuals.proportion,.$species.proportion, paired=TRUE))
w.summary <- proportions.test %>% tidy(w)

guild.means <- site.proportions.stat %>% group_by(trophic.guild) %>% summarise(mean.sp.prop = mean(species.proportion),
                                                                               sd.sp.prop = sd(species.proportion),
                                                                               mean.ind.prop = mean(individuals.proportion),
                                                                               sd.ind.prop = sd(individuals.proportion))
w.summary <- left_join(w.summary,guild.means)
w.summary <- w.summary[,c("trophic.guild","mean.sp.prop","sd.sp.prop","mean.ind.prop","sd.ind.prop","statistic","p.value")]

# 2 - dominance
# partial correlations, as in spencer2000

correlations <- dominance.data %>% group_by(trophic.guild) %>% summarize(correlation = cor(guild.richness,dominance))

################
# other: dominance-evenness?

# metrics data
all.metrics <- readr::read_delim(file = "./results/empirical_metrics.csv",delim = ";",col_types = list(col_character(),
                                                                                                          col_character(),
                                                                                                          col_double(),col_double(),col_double(),col_double(),
                                                                                                          col_double(),col_double(),col_double(),col_double(),
                                                                                                          col_double(),col_double(),col_double(),col_double()))

all.metrics$trophic.guild <- factor(all.metrics$trophic.guild, levels = c("plants","herbivores","omnivores","carnivores"))

#- guilds with more than X species, otherwise it does not make sense to compute SAD
minimum.richness <- 3
all.metrics <- subset(all.metrics, richness >= minimum.richness)
############
# natural or log-transformed abundances?
# empirical.metrics <- empirical.metrics[,c("site.ID","trophic.guild","hill.evenness.log","skewness.log")]
# names(empirical.metrics)[c(3,4)] <- c("hill.evenness","skewness")
############
empirical.metrics <- all.metrics[,c("site.ID","trophic.guild","hill.evenness","skewness")]
empirical.metrics <- droplevels(subset(empirical.metrics,trophic.guild != "plants"))
dominance.data <- dominance.data[,c("site.ID","trophic.guild","dominance","guild.richness")]

dominance.evenness <- left_join(dominance.data,empirical.metrics)

dom.ev.plot <- ggplot(dominance.evenness,aes(x = dominance, y = hill.evenness, fill = trophic.guild)) + 
  geom_point(aes(size = guild.richness), shape = 21) + 
  facet_grid(.~trophic.guild) +
  scale_color_manual(values = my.palette)+
  scale_fill_manual(values = my.palette)+
  
  xlab("dominance") + ylab("evenness") + 
  #DGC::theme_Publication() +
  theme(strip.background = element_blank())+#,strip.text.x = element_blank()) +
  guides(fill = FALSE) +
  NULL
dom.ev.plot
