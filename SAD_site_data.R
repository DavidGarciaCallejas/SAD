###############################
# clean up sites data
###############################

# read gentry sites data

gentry.sites <- readr::read_delim(file = "./data/abundance datasets/gentry/gentry_forest_transects_sites.csv",delim = ",")

# spatial and temporal resolution

gentry.sites$spatial.extent <- 1000 # m2, as stated in the documentation
gentry.sites$temporal.extent <- 0 # months

# habitat type

gentry.sites$habitat <- "MF/TF"

# clean the dataset

gentry.clean <- gentry.sites[,c("abbreviation","lat","lon","spatial.extent","temporal.extent","habitat")]
names(gentry.clean)[1] <- "site.ID"
names(gentry.clean)[2] <- "latitude"
names(gentry.clean)[3] <- "longitude"

gentry.clean$dataset <- "gentry"

##################################
##################################

# mammals data

mammals.sites <- readr::read_delim(file = "./data/abundance datasets/MCDB/MCDB_sites.csv",delim = ",")

# subset for abundance data

mammals.sites <- subset(mammals.sites,Abundance_data_present != "none" & Spatial_extent != "NULL" & Study_duration != "NULL")

# clean the dataset

mammals.clean <- mammals.sites[,c("Site_ID","Latitude","Longitude","Spatial_extent","Study_duration","Habitat_code")]
names(mammals.clean)[1] <- "site.ID"
names(mammals.clean)[2] <- "latitude"
names(mammals.clean)[3] <- "longitude"
names(mammals.clean)[4] <- "spatial.extent"
names(mammals.clean)[5] <- "temporal.extent"
names(mammals.clean)[6] <- "habitat"

mammals.clean$site.ID <- as.character(mammals.clean$site.ID)
mammals.clean$latitude <- as.numeric(mammals.clean$latitude)
mammals.clean$longitude <- as.numeric(mammals.clean$longitude)
mammals.clean$spatial.extent <- as.numeric(mammals.clean$spatial.extent)
mammals.clean$temporal.extent <- as.numeric(mammals.clean$temporal.extent)
mammals.clean$dataset <- "mammals"

####################
####################

all.sites <- bind_rows(gentry.clean,mammals.clean)

readr::write_delim(all.sites,"./data/abundance datasets/all_sites.csv",delim = ";")

