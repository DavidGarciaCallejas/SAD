
###############################
# clean up abundance data
###############################

# gentry abundances

gentry.abundances <- readr::read_delim(file = "./data/abundance datasets/gentry/gentry_forest_transects_counts.csv",delim = ",")

gentry.abundances <- gentry.abundances %>% group_by(site_code,species_id) %>% summarize(abundance = sum(count))
names(gentry.abundances)[1:2] <- c("site.ID","species.ID")
gentry.abundances$species.ID <- as.character(gentry.abundances$species.ID)
# mammals abundances

mammals.abundances <- readr::read_delim(file = "./data/abundance datasets/MCDB/MCDB_communities.csv",delim = ",")

mammals.abundances <- mammals.abundances[,c("Site_ID","Species_ID","Abundance")]
names(mammals.abundances) <- c("site.ID","species.ID","abundance")

mammals.abundances <- subset(mammals.abundances,abundance != "NULL")

mammals.abundances$site.ID <- as.character(mammals.abundances$site.ID)
mammals.abundances$abundance <- as.numeric(mammals.abundances$abundance)

##############
#wtf?
gentry.abundances <- as.data.frame(gentry.abundances)

all.abundances <- bind_rows(gentry.abundances,mammals.abundances)

readr::write_delim(all.abundances,"./data/abundance datasets/all_abundances.csv",delim = ";")



