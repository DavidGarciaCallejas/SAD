###############################
# clean up species data
###############################

# gentry species

gentry.sp <- readr::read_delim(file = "./data/abundance datasets/gentry/gentry_forest_transects_species.csv",delim = ",")

gentry.sp$trophic.guild <- "plants"
gentry.sp <- gentry.sp[,c("species_id","family","genus","species","trophic.guild")]
names(gentry.sp)[1] <- "species.ID"
gentry.sp$species.ID <- as.character(gentry.sp$species.ID)

# mammals species
# these are already classified in trophic guilds

mammals.sp <- readr::read_delim(file = "./results/mammals_trophic_classification.csv",delim = ";")

all.species <- bind_rows(gentry.sp,mammals.sp)

# write species
readr::write_delim(all.species,"./data/abundance datasets/all_species.csv",delim = ";")


