#####################
# group trophic positions of a food web niche model into trophic levels
#####################

# model.data <- readr::read_delim(file = "./results/abundances_niche_model_DD.csv",delim = ";")
model.data <- readr::read_delim(file = "./results/abundances_niche_model_DP.csv",delim = ";")
# model.data <- readr::read_delim(file = "./results/abundances_niche_model_RF.csv",delim = ";")

richness.levels <- unique(model.data$richness)
connectance.levels <- unique(model.data$connectance)
apportionment.levels <- unique(model.data$niche.apport)
replicates <- max(model.data$replicate)

model.data <- model.data[complete.cases(model.data),]
model.data$trophic.guild <- "-"

# i.richness <- 1
# i.connectance <- 1
# i.apport <- 1
# i.replicate <- 1

for(i.richness in 1:length(richness.levels)){
  for(i.connectance in 1:length(connectance.levels)){
    for(i.apport in 1:length(apportionment.levels)){
      for(i.replicate in 1:replicates){
        
        my.index <- which(model.data$richness == richness.levels[i.richness] & 
                                model.data$connectance == connectance.levels[i.connectance] &
                                model.data$niche.apport == apportionment.levels[i.apport] &
                                model.data$replicate == i.replicate)
        
        for(i.sp in 1:length(my.index)){
          my.trophic.position <- model.data$trophic.position[my.index[i.sp]]
          if(my.trophic.position == 1){
            model.data$trophic.guild[my.index[i.sp]] <- "primary Producers"
          }else if(my.trophic.position > 1 & my.trophic.position <= 2.3){
            model.data$trophic.guild[my.index[i.sp]] <- "herbivores"
          }else if(model.data$plant.consumer[my.index[i.sp]] == TRUE){
            model.data$trophic.guild[my.index[i.sp]] <- "omnivores"
          }else{
            model.data$trophic.guild[my.index[i.sp]] <- "carnivores"
          }
        }
        
        #Williams and Martinez 2008 have these categories
        # my.data$trophic.guild <- "Omnivores"
        # my.data$trophic.guild[my.data$trophic.position == 1] <- "Primary Producers"
        # my.data$trophic.guild[my.data$trophic.position >= 2 & my.data$trophic.position < 2.3] <- "Herbivores"
        # my.data$trophic.guild[my.data$trophic.position > 2.7 & my.data$trophic.position < 3.3] <- "Carnivores"
        # my.data$trophic.guild[my.data$trophic.position > 3.7 & my.data$trophic.position < 4.3] <- "Carnivores"
        # my.data$trophic.guild[my.data$trophic.position > 4.7 & my.data$trophic.position < 5.3] <- "Carnivores"
        # my.data$trophic.guild[my.data$trophic.position > 5.7 & my.data$trophic.position < 6.3] <- "Carnivores"
        # my.data$trophic.guild[my.data$trophic.position > 6.7 & my.data$trophic.position < 7.3] <- "Carnivores"
        # my.data$trophic.guild[my.data$trophic.position > 7.7 & my.data$trophic.position < 8.3] <- "Carnivores"
        
        # model.data$trophic.guild[my.index] <- "Omnivores"
        # my.data$trophic.guild[my.data$trophic.position == 1] <- "Primary Producers"
        # my.data$trophic.guild[my.data$trophic.position >= 2 & my.data$trophic.position <= 2.3] <- "Herbivores"
        # my.data$trophic.guild[my.data$trophic.position > 2.3 & my.data$plant.consumer == FALSE] <- "Carnivores"
        
      }# for i.replicate
    }# for i.apport
  }# for i.connectance
}# for i.richness

model.data <- model.data[,c("richness","connectance","niche.apport","replicate","ID","M","N_predicted","trophic.position","trophic.guild")]
# readr::write_delim(x = model.data,path = "./results/abundances_niche_model_trophic_guild_DD.csv",delim = ";")
readr::write_delim(x = model.data,path = "./results/abundances_niche_model_trophic_guild_DP.csv",delim = ";")
# readr::write_delim(x = model.data,path = "./results/abundances_niche_model_trophic_guild_RF.csv",delim = ";")

