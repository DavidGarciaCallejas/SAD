###############################
# mammals classification in trophic guilds
# this version: herbivores, omnivores, carnivores
###############################

mammals.guilds <- readr::read_delim(file = "./data/trophic_classification/mammals.csv",delim = ";")

trophic.indexes <- mammals.guilds[,1:3]
names(trophic.indexes) <- c("trophic.ID","full species","family")
trophic.indexes$genus <- character(nrow(trophic.indexes))
trophic.indexes$species <- character(nrow(trophic.indexes))
trophic.genus <- strsplit(trophic.indexes$`full species`," ")
for(i.sp in 1:nrow(trophic.indexes)){
  trophic.indexes$genus[i.sp] <- trophic.genus[[i.sp]][1]
  trophic.indexes$species[i.sp] <- trophic.genus[[i.sp]][2]
}
trophic.indexes <- trophic.indexes[,c("trophic.ID","family","genus","species")]

mammals.guilds <- mammals.guilds[,c("MSW3_ID",
                                    "Diet-Inv",
                                    "Diet-Vend",
                                    "Diet-Vect",
                                    "Diet-Vfish",
                                    "Diet-Vunk",
                                    "Diet-Scav",
                                    "Diet-Fruit",
                                    "Diet-Nect",
                                    "Diet-Seed",
                                    "Diet-PlantO",
                                    "BodyMass-Value")]

#####
mammals.guilds$plants <- mammals.guilds$`Diet-PlantO` + 
  mammals.guilds$`Diet-Seed` + 
  mammals.guilds$`Diet-Fruit` + 
  mammals.guilds$`Diet-Nect`
mammals.guilds$carn <- mammals.guilds$`Diet-Vend` + 
  mammals.guilds$`Diet-Vect` + 
  mammals.guilds$`Diet-Vfish` + 
  mammals.guilds$`Diet-Vunk` + 
  mammals.guilds$`Diet-Scav` + 
  mammals.guilds$`Diet-Inv`

mammals.guilds$trophic.guild <- character(nrow(mammals.guilds))

for(i.sp in 1:nrow(mammals.guilds)){
  if(mammals.guilds$plants[i.sp] > 70){
    mammals.guilds$trophic.guild[i.sp] <- "herbivores"
  }else if(mammals.guilds$carn[i.sp] > 70){
    mammals.guilds$trophic.guild[i.sp] <- "carnivores"
  }else{
    mammals.guilds$trophic.guild[i.sp] <- "omnivores"
  }
}# for each sp

############
# same for genus, to extrapolate unclassified species

# group by genus
mammals.guilds$genus <- character(nrow(mammals.guilds))
mammals.guilds$species <- character(nrow(mammals.guilds))
for(i.sp in 1:nrow(mammals.guilds)){
  mammals.guilds$genus[i.sp] <- trophic.genus[[i.sp]][1]
  mammals.guilds$species[i.sp] <- trophic.genus[[i.sp]][2]
}
mammals.guilds.tidy <- mammals.guilds %>% gather(key = diet,value = perc,`Diet-Inv`:`Diet-PlantO`)
mammals.guilds.tidy <- mammals.guilds.tidy[,c("MSW3_ID","genus","species","diet","perc")]
genus.diet <- mammals.guilds.tidy %>% group_by(genus,diet) %>% summarise(mean.perc = mean(perc))
genus.diet <- spread(genus.diet,diet,mean.perc)

#####
genus.diet$plants <- genus.diet$`Diet-PlantO` + 
  genus.diet$`Diet-Seed` + 
  genus.diet$`Diet-Fruit` + 
  genus.diet$`Diet-Nect`
genus.diet$carn <- genus.diet$`Diet-Vend` + 
  genus.diet$`Diet-Vect` + 
  genus.diet$`Diet-Vfish` + 
  genus.diet$`Diet-Vunk` + 
  genus.diet$`Diet-Scav` + 
  genus.diet$`Diet-Inv`

genus.diet$trophic.guild <- character(nrow(genus.diet))

for(i.sp in 1:nrow(genus.diet)){
  if(genus.diet$plants[i.sp] > 70){
    genus.diet$trophic.guild[i.sp] <- "herbivores"
  }else if(genus.diet$carn[i.sp] > 70){
    genus.diet$trophic.guild[i.sp] <- "carnivores"
  }else{
    genus.diet$trophic.guild[i.sp] <- "omnivores"
  }
}# for each sp

###
# 2 - assign the trophic guilds to the mammal species of the SAD datasets

mammals.sp <- readr::read_delim(file = "./data/abundance datasets/MCDB/MCDB_species.csv",delim = ",")
mammals.sp <- mammals.sp[,c("Species_ID","Family","Genus","Species")]
names(mammals.sp) <- c("species.ID","family","genus","species")

mammals.guild.keys <- mammals.guilds[,c("MSW3_ID","genus","species","trophic.guild")]

mammals.sp <- left_join(mammals.sp,mammals.guild.keys)

# unclassified ones
unclassified <- subset(mammals.sp,is.na(trophic.guild))
unclassified$genus <- gsub(" ", "", unclassified$genus)

unclassified$trophic.guild <- genus.diet$trophic.guild[match(unclassified$genus,genus.diet$genus)]

mammals.sp$trophic.guild[is.na(mammals.sp$trophic.guild)] <- unclassified$trophic.guild[match(mammals.sp$species.ID[is.na(mammals.sp$trophic.guild)],
                                                                                              unclassified$species.ID)]
# still two unclassified: soricidae and tayassuidae
# classify manually (beware the clusters made by the algorithm)
# or discard
mammals.sp <- subset(mammals.sp,!is.na(trophic.guild))
#########

mammals.sp <- mammals.sp[,c("species.ID","family","genus","species","trophic.guild")]
readr::write_delim(mammals.sp,path = "./results/mammals_trophic_classification.csv",delim = ";")
