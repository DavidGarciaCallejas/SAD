
#####################
# goodness of fit of the tuesday lake data
#####################

# alpha value: do I want to use M^1 or M^alpha??
alpha.metabolism <- 1
#####################

CohenAbundance <- function(species.mass,prey.data,type = "N", alpha = 1){

    if(nrow(prey.data)>0){
      sum.biomass <- sum((prey.data$N*(prey.data$M^alpha))/prey.data$V)
      if(type == "N"){
        predicted.abundance <- sum.biomass/species.mass
      }else{
        predicted.abundance <- sum.biomass
      }
    }else{
      stop("CohenAbundance: error in provided data, consumer species found without prey species")
    }
  predicted.abundance
}# function

#####################
#####################
Vulnerability <- function(interaction.matrix,species.data,return.type = "net"){
  species.consumers <- numeric(nrow(species.data))
  species.potential.consumers <- numeric(nrow(species.data))
  for(i.sp in 1:nrow(species.data)){
    # number of potential consumers
    species.potential.consumers[i.sp] <- sum(species.data$trophic.level == (species.data$trophic.level[i.sp]+1))
    # number of consumers
    species.consumers[i.sp] <- sum(interaction.matrix[i.sp,] == -1)
  }# for i.sp
  
  if(return.type == "net"){
    return(species.consumers)
  }else if(return.type == "ratio"){
    return(species.consumers/species.potential.consumers)
  }# if-else
  
}

#####################
#####################
# read empirical datasets and tidy them to
# interaction matrix + data frame

# TUESDAY LAKE

# 1 - read species data
tuesday.lake.data.1984 <- read.table(file = "./data/abundance datasets/tuesday lake/Tuesday_lake_1984_sp.csv",header = T,sep = ";",stringsAsFactors = F)
tuesday.lake.data.1986 <- read.table(file = "./data/abundance datasets/tuesday lake/Tuesday_lake_1986_sp.csv",header = T,sep = ";",stringsAsFactors = F)

# in the original names are the actual units
names(tuesday.lake.data.1984)[8] <- "biomass(kg)"
names(tuesday.lake.data.1984)[9] <- "abundance"

names(tuesday.lake.data.1986)[8] <- "biomass(kg)"
names(tuesday.lake.data.1986)[9] <- "abundance"

# clean up dataframe
tuesday.lake.data.1984 <- tuesday.lake.data.1984[,c("ID","Genus.species","abundance","biomass(kg)","trophic.position.grouping")]
tuesday.lake.data.1986 <- tuesday.lake.data.1986[,c("ID","Name","abundance","biomass(kg)","trophic.position.grouping")]

names(tuesday.lake.data.1984) <- c("ID","name","N","M","trophic.level")
names(tuesday.lake.data.1986) <- c("ID","name","N","M","trophic.level")

tuesday.lake.data.1984$trophic.level[tuesday.lake.data.1984$trophic.level == "basal"] <- 1
tuesday.lake.data.1984$trophic.level[tuesday.lake.data.1984$trophic.level == "intermediate"] <- 2
tuesday.lake.data.1984$trophic.level[tuesday.lake.data.1984$trophic.level == "top"] <- 3
tuesday.lake.data.1984$trophic.level <- as.numeric(tuesday.lake.data.1984$trophic.level)

tuesday.lake.data.1986$trophic.level[tuesday.lake.data.1986$trophic.level == "basal"] <- 1
tuesday.lake.data.1986$trophic.level[tuesday.lake.data.1986$trophic.level == "intermediate"] <- 2
tuesday.lake.data.1986$trophic.level[tuesday.lake.data.1986$trophic.level == "top"] <- 3
tuesday.lake.data.1986$trophic.level <- as.numeric(tuesday.lake.data.1986$trophic.level)

# tuesday.lake.data.1984 <- subset(tuesday.lake.data.1984,trophic.level != "-1")
# tuesday.lake.data.1986 <- subset(tuesday.lake.data.1986,trophic.level != "-1")

# read links data and generate interaction matrix
links.1984 <- read.table(file = "./data/abundance datasets/tuesday lake/tuesday_lake_1984_links.csv",header = T,sep = ";")
links.1986 <- read.table(file = "./data/abundance datasets/tuesday lake/tuesday_lake_1986_links.csv",header = T,sep = ";")

tuesday.lake.matrix.1984 <- matrix(0,nrow(tuesday.lake.data.1984),nrow(tuesday.lake.data.1984))
for(i.link in 1:nrow(links.1984)){
  
  tuesday.lake.matrix.1984[links.1984$Resource.ID[i.link],links.1984$Consumer.ID[i.link]] <- -1
  tuesday.lake.matrix.1984[links.1984$Consumer.ID[i.link],links.1984$Resource.ID[i.link]] <- 1
  
}# for i.link

tuesday.lake.matrix.1986 <- matrix(0,nrow(tuesday.lake.data.1986),nrow(tuesday.lake.data.1986))
for(i.link in 1:nrow(links.1986)){
  
  tuesday.lake.matrix.1986[links.1986$Resource.ID[i.link],links.1986$Consumer.ID[i.link]] <- -1
  tuesday.lake.matrix.1986[links.1986$Consumer.ID[i.link],links.1986$Resource.ID[i.link]] <- 1
  
}# for i.link

#####################
#####################
# calculate the predicted abundances for each dataset

datasets <- c("Tuesday Lake - 1984","Tuesday Lake - 1986")

# observed-predicted
obs.pred.data <- data.frame(species = integer(1e5),
                            dataset = character(1e5),
                            trophic.level = numeric(1e5),
                            V = integer(1e5),M = numeric(1e5),
                            B_obs = numeric(1e5),B_pred_single = numeric(1e5),B_pred_full = numeric(1e5),
                            N_obs = numeric(1e5),N_pred_single = numeric(1e5),N_pred_full = numeric(1e5),
                            stringsAsFactors = F)
count <- 1

for(i.dataset in 1:length(datasets)){
  
  if(datasets[i.dataset] == "Tuesday Lake - 1984"){
    orig.matrix <- tuesday.lake.matrix.1984
    interaction.matrix <- matrix(data = 0,nrow = nrow(tuesday.lake.matrix.1984),ncol = nrow(tuesday.lake.matrix.1984))
    my.species.data <- tuesday.lake.data.1984
  }else if(datasets[i.dataset] == "Tuesday Lake - 1986"){
    orig.matrix <- tuesday.lake.matrix.1986
    interaction.matrix <- matrix(data = 0,nrow = nrow(tuesday.lake.matrix.1986),ncol = nrow(tuesday.lake.matrix.1986))
    my.species.data <- tuesday.lake.data.1986
  }
  
  #################
  #################
  
      # first trophic level, predicted = observed 
    for(i.sp in 1:nrow(my.species.data)){
      
      
      if(my.species.data$trophic.level[i.sp] %in% c(-1,1)){
        obs.pred.data$species[count] <- my.species.data$ID[i.sp]
        obs.pred.data$dataset[count] <- datasets[i.dataset]
        obs.pred.data$trophic.level[count] <- my.species.data$trophic.level[i.sp]
        obs.pred.data$V[count] <- sum(orig.matrix[i.sp,] == -1)
        obs.pred.data$M[count] <- my.species.data$M[i.sp]
        obs.pred.data$N_obs[count] <- my.species.data$N[i.sp]
        obs.pred.data$B_obs[count] <- my.species.data$M[i.sp] * my.species.data$N[i.sp]
        
        # if first trophic level, no predicted values...
        obs.pred.data$N_pred_single[count] <- obs.pred.data$N_obs[count]
        obs.pred.data$B_pred_single[count] <- obs.pred.data$B_obs[count]
        
        obs.pred.data$N_pred_full[count] <- obs.pred.data$N_obs[count]
        obs.pred.data$B_pred_full[count] <- obs.pred.data$B_obs[count]
        
        count <- count + 1
      }
    }
    
    ######################
    # upper trophic levels
    
    left.sp <- my.species.data$ID[my.species.data$trophic.level > 1]
    sp.index <- 1
    
    # which species can we calculate
    while(length(left.sp)>0 & !(sp.index > length(left.sp))){
      my.sp <- left.sp[sp.index]
      prey.sp <- which(orig.matrix[,my.sp] == -1)
      # any sp with all prey sp calculated?
      while(sum(prey.sp %in% left.sp)>0 & !(sp.index > length(left.sp))){
        sp.index <- sp.index + 1
        my.sp <- left.sp[sp.index]
        prey.sp <- which(orig.matrix[,my.sp] == -1)
      }
      # if any, calculate it
      if(sum(prey.sp %in% left.sp) == 0){
        my.index <- which(my.species.data$ID == my.sp)
        obs.pred.data$species[count] <- my.sp
        obs.pred.data$dataset[count] <- datasets[i.dataset]
        obs.pred.data$trophic.level[count] <- my.species.data$trophic.level[my.index]
        obs.pred.data$V[count] <- sum(orig.matrix[my.sp,] == -1)
        obs.pred.data$M[count] <- my.species.data$M[my.index]
        obs.pred.data$N_obs[count] <- my.species.data$N[my.index]
        obs.pred.data$B_obs[count] <- my.species.data$M[my.index] * my.species.data$N[my.index]
        
        # predicted prey abundances
        prey.data <- data.frame(ID = prey.sp,
                                M = obs.pred.data$M[obs.pred.data$species %in% prey.sp & obs.pred.data$dataset == datasets[i.dataset]],
                                N = obs.pred.data$N_pred_full[obs.pred.data$species %in% prey.sp & obs.pred.data$dataset == datasets[i.dataset]],
                                V = 0)
        # vulnerability of the prey species
        for(i.prey in 1:length(prey.sp)){
          prey.data$V[i.prey] <- sum(orig.matrix[prey.sp[i.prey],] == -1)
        }
        
        # observed prey abundances
        prey.observed.data <- my.species.data[my.species.data$ID %in% prey.sp,]
        prey.observed.data$V <- prey.data$V[match(prey.observed.data$ID,prey.data$ID)]
        
        # predicted abundance using observed prey abundances
        obs.pred.data$N_pred_single[count] <- CohenAbundance(species.mass = obs.pred.data$M[count], prey.data = prey.observed.data, type = "N", alpha = alpha.metabolism)
        obs.pred.data$B_pred_single[count] <- CohenAbundance(species.mass = obs.pred.data$M[count], prey.data = prey.observed.data, type = "B", alpha = alpha.metabolism)
        
        # predicted abundance using predicted prey abundances
        obs.pred.data$N_pred_full[count] <- CohenAbundance(species.mass = obs.pred.data$M[count], prey.data = prey.data, type = "N", alpha = alpha.metabolism)
        obs.pred.data$B_pred_full[count] <- CohenAbundance(species.mass = obs.pred.data$M[count], prey.data = prey.data, type = "B", alpha = alpha.metabolism)
        
        # update results counter
        count <- count + 1
        # remove this sp from the list of left species
        left.sp <- left.sp[-sp.index]
        # restart sp.index counter
        sp.index <- 1
      }# if valid sp
    }# while
}# for each dataset

obs.pred.data <- subset(obs.pred.data,species != 0)

if(alpha.metabolism != 1){
  my.file.name <- "./results/cohen_tuesday_lake_abundances_with_alpha.csv"
}else{
  my.file.name <- "./results/cohen_tuesday_lake_abundances.csv"
}

readr::write_delim(x = obs.pred.data,path = my.file.name,delim = ";")

