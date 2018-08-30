
############################
# calculate variance, skewness and Hill evenness for abundance datasets
############################

# abundances data
abundances <- readr::read_delim(file = "./data/abundance datasets/all_abundances.csv",delim = ";",col_types = list(col_character(),col_character(),col_double()))
abundances <- subset(abundances, !is.na(abundance))

# species data
species.data <- readr::read_delim(file = "./data/abundance datasets/all_species.csv",
                                  delim = ";",
                                  col_types = list(col_character(),
                                                   col_character(),
                                                   col_character(),
                                                   col_character(),col_character()))

############
# calculate indices for every trophic guild of each site
abundances$trophic.guild <- species.data$trophic.guild[match(abundances$species.ID,species.data$species.ID)]

# some cleaning
abundances <- subset(abundances,!is.na(trophic.guild))
abundances <- subset(abundances, abundance > 0)

# transform data for consistency
# abundances < 1 -> 1.001
# this way log.data can be computed and hill.diversity function does not raise errors
# note that in empirical datasets there should not be fractional abundances, 
# but abundances = 1 are certain to occur.

abundances$abundance[abundances$abundance <= 1] <- 1.001
abundances$log.abundance <- log(abundances$abundance,2)

site.metrics <- abundances %>% group_by(site.ID,trophic.guild) %>% summarise(variance = var(abundance),
                                                                             mad = stats::mad(abundance, constant = 1),
                                                                             skewness = mc(abundance),
                                                                             hill.diversity = hill.diversity(abundance),
                                                                             richness = sum(abundance>0))
site.metrics$hill.evenness <- site.metrics$hill.diversity/site.metrics$richness

### now, with logs
site.metrics.log <- abundances %>% group_by(site.ID,trophic.guild) %>% summarise(variance.log = var(log.abundance),
                                                                                 mad.log = stats::mad(log.abundance, constant = 1),
                                                                                 skewness.log = mc(log.abundance),
                                                                                 hill.diversity.log = hill.diversity(log.abundance),
                                                                                 richness.log = n())
site.metrics.log$hill.evenness.log <- site.metrics.log$hill.diversity.log/site.metrics$richness

################################
# a first test
# ev.plot <- ggplot(site.metrics,aes(x = trophic.guild,y = hill.evenness, group = trophic.guild)) +
#   geom_boxplot()
# ev.plot
################################

site.metrics <- left_join(site.metrics,site.metrics.log)

readr::write_delim(x = site.metrics,path = "./results/empirical_metrics.csv",delim = ";",append = F)
