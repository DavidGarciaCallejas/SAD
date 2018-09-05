# theoretical distributions of abundances in increasing trophic levels
# given body sizes of all levels, and abundances of the basal level

#####################
# auxiliary functions

# niche model version from Kleka 2014 plos one
NicheModelv2<-function(S,C,mass.range){ # S = the number of species, C = target connectance, mass.range = the range of body masses in the food web (w_max/w_min)
  
  C.web=2 # C.web = actual connectance, set C.web to a nonsensical value to initiate the while loop
  
  while((C.web>1.03*C)|(C.web<0.97*C)){ # C.web should be +- 3% of the target value of C
    
    con<-rep(0,times=S) # vector containing 0 for disconnected species
    
    while(length(con[con==0])>0){ # while some elements of con are zero, generate a new web until there are no disconnected species
      
      beta<-(1-2*C)/(2*C)
      species.niche<-sort(runif(n=S,min=0,max=1)) # generate niche values
      
      body.masses<-10^(mass.range*species.niche) # calculate body masses from niche values
      
      species.ri<-c(1-(1-runif(S))^(1/beta))*species.niche
      species.ri[1]<-0 #the smallest species is always a primary producer, it has zero niche width
      
      for (n in 1:S){
        while (species.ri[n]/2>species.niche[n]) {species.ri[n]<-(1-(1-runif(1))^(1/beta))} # species r_i has to satisfy this condition
      }
      
      species.centre<-sapply(1:S,function (x) runif(1,min=species.ri[x]/2,max=species.niche[x])) # the centre of the feeding niche
      
      web<-matrix(0,nrow=S,ncol=S) # initiate the food web matrix
      
      for (j in 1:S){
        for (i in 1:S){
          if ((species.niche[i]<(species.centre[j]+species.ri[j]/2)) & (species.niche[i]>(species.centre[j]-species.ri[j]/2))) {web[i,j]<-1}
          # fill in the food web matrix
        }
      }
      
      web.zero.diagonal<-web*(matrix(1,nrow=S,ncol=S)-diag(S)) # set all values on the diagonal to zero, so that purely cannibalistic species are counted as disconnected
      con<-colSums(web.zero.diagonal)+rowSums(web.zero.diagonal) # disconnected species have a value of con=0
    }
    C.web=sum(web)/S^2
  }
  
  return(list("model"=c("original niche model"),"body.masses"=body.masses,"web"=web,"S"=S,"C"=C,"mass.range"=mass.range))
}

#####################
# the niche model, own version
NicheModel <- function(S,connectance,niche = NULL){
  
  interaction.matrix <- matrix(0,nrow=S,ncol=S)
  
  # niche axis
  if(is.null(niche)){
    my.niche <- sort(runif(S,0,1))
  }else{
    my.niche <- niche
  }  
  # Range
  beta = 1/(2*connectance)-1
  range = rbeta(S,1,beta)*my.niche
  
  #Centroid
  centroid = runif(S,range/2,my.niche)
  
  # Make the first niche a producer
  range[which.min(my.niche)] = .0000001
  
  # Evaluate the matrix of links
  
  low = centroid - range/2
  high = centroid + range/2
  
  for(s1 in 1:S){
    for(s2 in 1:S){
      if(s1 != s2){
        if(low[s1] < my.niche[s2] && high[s1] > my.niche[s2]){ 
          interaction.matrix[s2,s1] = 1    	  
        }# if
      }# if
    }# for
  }# for
  if(is.null(niche)){
    return(list(interaction.matrix,my.niche))
  }else{
    return(interaction.matrix)  
  } 
}

#####################
# network parameters
connectance.levels <- c(0.1,0.2,0.5)
richness.levels <- c(50,100,200)

# basal abundances distribution
shape.weibull.min <- 0.15
shape.weibull.max <- 0.2
scale.weibull <- 4.7 * 100

# other parameters
# Lindeman's ten percent law
# trebilco et al. 2013
conversion.efficiency.min <- 0.1
conversion.efficiency.max <- 0.15

# simulations
replicates <- 1000

#######################
sim.results <- NULL

# 
# i.richness <- 1
# i.connectance <- 1
# i.replicate <- 1

for(i.richness in 1:length(richness.levels)){
  for(i.connectance in 1:length(connectance.levels)){
    for(i.replicate in 1:replicates){
      
      # results dataframe
      sp.data <- data.frame(richness = richness.levels[i.richness],
                            connectance = connectance.levels[i.connectance],
                            replicate = i.replicate,
                            ID = 1:richness.levels[i.richness],
                            trophic.position = 0,
                            trophic.guild = "-",
                            M = 0,
                            N_predicted = 0,
                            # auxiliary
                            complete = FALSE)
      
      # niche axis
      my.niche <- sort(runif(richness.levels[i.richness],0,1))
      # generate network with the niche model
      interaction.matrix <- NicheModel(richness.levels[i.richness],connectance.levels[i.connectance],my.niche)
      
      # body sizes - check
      sp.data$M <- my.niche
      
      # do not allow closed loops, for example sp1 feeding on sp2 and sp2 feeding on sp1
      # randomly remove one of the links
      # longer loops are removed below
      for(i.sp in 1:nrow(sp.data)){
        for(j.sp in 1:nrow(sp.data)){
          if(interaction.matrix[i.sp,j.sp] == 1 & interaction.matrix[j.sp,i.sp] == 1){
            to.remove <- sample(1:2,1)
            if(to.remove == 1){
              interaction.matrix[i.sp,j.sp] <- 0
            }else{
              interaction.matrix[j.sp,i.sp] <- 0
            }
          }
        }
      }
      
      basal <- which(colSums(interaction.matrix) == 0)
      # check this part
      sp.data$N_predicted[basal] <- rweibull(length(basal),
                                              shape = runif(1,shape.weibull.min,shape.weibull.max),
                                              scale = scale.weibull)
      # for now!
      # sp.data$N_predicted[sp.data$N_predicted < 1] <- 1
      
      ##################
      # generate predicted abundances
      # first, propagate from basal species
      for(i.basal in 1:length(basal)){
        
        # if the basal sp has predators
        if(sum(interaction.matrix[basal[i.basal],])>0){
          
          resource.available <- sp.data$M[basal[i.basal]] * sp.data$N_predicted[basal[i.basal]] * runif(1,conversion.efficiency.min,conversion.efficiency.max)
          predator.set <- data.frame(predator.id = which(interaction.matrix[basal[i.basal],] == 1), specialization = 0)
          if(nrow(predator.set)>1){
            # resource is partitioned according to the degree of specialization
            # and this is given by the number of links of the predators
            for(i.predator in 1:nrow(predator.set)){
              predator.set$specialization[i.predator] <- sum(interaction.matrix[,predator.set$predator.id[i.predator]] == 1)
            }# for each predator
            
            # the more specialized predators are the best competitors
            predator.set <- arrange(predator.set,specialization)
            equal.specialization <- predator.set$predator.id[which(predator.set$specialization %in% predator.set$specialization[which(duplicated(predator.set$specialization))])]
            predator.order <- predator.set$predator.id
            if(length(equal.specialization)>0){
              equal.order <- which(predator.order %in% equal.specialization)
              new.order <- sample(equal.order)
              predator.order[equal.order] <- predator.order[new.order]
            }
            
            # dominance preemption, sugihara 1980, tokeshi 1990
            resource.partition <- nicheApport::dominancePreemp(resource.available,nrow(predator.set))
            # if i ever use dominance decay, beware that it fails for n=2
            # predator.set <- arrange(predator.set,specialization)
            # predator.set$assigned <- resource.partition
            
            sp.data$N_predicted[sp.data$ID %in% predator.order] <- sp.data$N_predicted[sp.data$ID %in% predator.order] + 
              resource.partition
          }else if(nrow(predator.set) == 1){
            sp.data$N_predicted[sp.data$ID == predator.set] <- resource.available
          }#else if(nrow(predator.set) == 0){
          # no predators, so nothing to do
          #}
        }# if it has predators  
        
        sp.data$complete[basal[i.basal]] <- TRUE
      }# for i.basal
      
      # second, loop through the rest of species, propagating the biomass flows
      # and substituting interactions if a loop is found
      while(!all(sp.data$complete)){
        
        current.complete <- sum(sp.data$complete)
        # after basal species, all species that only feed on basal ones should be complete
        for(i.sp in 1:nrow(sp.data)){
          
          if(sp.data$complete[i.sp] == FALSE){
            
            prey.set <- which(interaction.matrix[,i.sp] == 1)
            if(length(prey.set)>0){
              complete.sp.match <- match(prey.set,sp.data$ID[sp.data$complete == TRUE])
              # if all preys are complete, 
              # this sp is also complete and can "propagate" up in the network
              if(!anyNA(complete.sp.match)){
                sp.data$complete[i.sp] <- TRUE
                remaining <- sum(sp.data$complete == F)
                # print(paste("sp:",i.sp," complete, ", remaining, " remaining", sep=""))
                # propagate its abundance, but only if it has predators
                if(sum(interaction.matrix[i.sp,])>0){
                  
                  resource.available <- sp.data$M[i.sp] * sp.data$N_predicted[i.sp] * runif(1,conversion.efficiency.min,conversion.efficiency.max)
                  predator.set <- data.frame(predator.id = which(interaction.matrix[i.sp,] == 1), specialization = 0)
                  if(nrow(predator.set)>1){
                    # resource is partitioned according to the degree of specialization
                    # and this is given by the number of links of the predators
                    for(i.predator in 1:nrow(predator.set)){
                      predator.set$specialization[i.predator] <- sum(interaction.matrix[,predator.set$predator.id[i.predator]] == 1)
                    }# for each predator
                    
                    # dominance preemption, sugihara 1980, tokeshi 1990
                    resource.partition <- nicheApport::dominancePreemp(resource.available,nrow(predator.set))
                    # if i ever use dominance decay, beware that it fails for n=2
                    # predator.set <- arrange(predator.set,specialization)
                    # predator.set$assigned <- resource.partition
                    
                    sp.data$N_predicted[sp.data$ID %in% predator.set$predator.id] <- sp.data$N_predicted[sp.data$ID %in% predator.set$predator.id] + 
                      resource.partition
                  }else if(nrow(predator.set) == 1){
                    sp.data$N_predicted[sp.data$ID == predator.set$predator.id] <- sp.data$N_predicted[sp.data$ID == predator.set$predator.id] + 
                      resource.available
                  }
                }# does it have predators?
              }# if any NA
            }# if predator
            
          } # if not complete
        }# for each sp
        after.complete <- sum(sp.data$complete)
        
        # are there loops without primary producers?
        # this happens if not all species are "complete"
        # in two succesive rounds
        # if this is the case, try to "break" the loop
        # by removing an interaction and reassigning it to a primary producer-consumer interaction
        if(after.complete == current.complete){
          # print(paste("entering reassigning phase..."))
          incomplete.sp <- sp.data[sp.data$complete == FALSE,]
          
          initial.sp <- incomplete.sp$ID[1]  
          my.prey <- which(interaction.matrix[,initial.sp] == 1)
          incomplete.prey <- sp.data$ID[sp.data$ID %in% my.prey & sp.data$complete == FALSE]
          loop.sp <- c(initial.sp)
          
          while(anyDuplicated(loop.sp) == 0){
            loop.sp <- c(loop.sp,incomplete.prey[1])
            my.sp <- incomplete.prey[1]
            my.prey <- which(interaction.matrix[,my.sp] == 1)
            incomplete.prey <- sp.data$ID[sp.data$ID %in% my.prey & sp.data$complete == FALSE]
          }
          
          selected.consumer <- sample.vec(loop.sp,1)#loop.sp[anyDuplicated(loop.sp)]
          
          # in order to increase artificially the number of herbivores, 
          # substitute all the links of the selected sp for basal-herbivore links
          # this may alter the structural patterns of the niche model, check below
          
          num.links <- sum(interaction.matrix[,selected.consumer])
          # use the vector "basal" from a few lines above
          interaction.matrix[,selected.consumer] <- 0
          num.links <- ifelse(num.links>length(basal),length(basal),num.links)
          
          new.resources <- sample.vec(basal,num.links)
          
          interaction.matrix[new.resources,selected.consumer] <- 1
          
          # alternatively, substitute only one link
          # selected.resource <- which(interaction.matrix[,selected.consumer] == 1)
          # selected.resource <- sp.data$ID[sp.data$ID %in% selected.resource & sp.data$complete == FALSE]
          # selected.resource <- sample.vec(selected.resource,1)
          # 
          # interaction.matrix[selected.resource,selected.consumer] <- 0
          # 
          # # substitute the removed link with a new link from a basal resource
          # new.resource <- which(colSums(interaction.matrix) == 0)
          # new.resource <- new.resource[which(interaction.matrix[new.resource,selected.consumer] == 0)]
          # new.resource <- sample.vec(new.resource,1)
          # 
          # interaction.matrix[new.resource,selected.consumer] <- 1
        }
        
      }# while species not completed
      
      ###########################
      # trophic level assignment
      cheddar.dataframe <- sp.data
      names(cheddar.dataframe)[which(names(cheddar.dataframe) == "ID")] <- "node"
      cheddar.dataframe <- cheddar.dataframe[,c("node","M")]
      cheddar.dataframe$node <- as.character(cheddar.dataframe$node)
      cheddar.dataframe$node <- paste("sp",cheddar.dataframe$node,sep="")
      cheddar.links <- expand.grid(1:nrow(sp.data),1:nrow(sp.data))
      names(cheddar.links)[1] <- "resource"
      names(cheddar.links)[2] <- "consumer"
      cheddar.links$maintain <- TRUE
      for(i.link in 1:nrow(cheddar.links)){
        if(interaction.matrix[cheddar.links$resource[i.link],cheddar.links$consumer[i.link]] == 0){
          cheddar.links$maintain[i.link] <- FALSE
        }
      }
      cheddar.links <- subset(cheddar.links,maintain == TRUE)
      cheddar.links <- cheddar.links[,1:2]
      cheddar.links$resource <- paste("sp",cheddar.links$resource,sep="")
      cheddar.links$consumer <- paste("sp",cheddar.links$consumer,sep="")
      
      cheddar.community <- Community(cheddar.dataframe,properties = list(title="foo",M.units="arbitrary"),trophic.links = cheddar.links)
      trophic.levels.community <- PreyAveragedTrophicLevel(cheddar.community,include.isolated = F)
      sp.data$trophic.position <- trophic.levels.community
      
      # hist(sp.data$trophic.position,breaks = 100)
      
      #Williams and Martinez 2008 have these categories
      sp.data$trophic.guild <- "Omnivores"
      sp.data$trophic.guild[sp.data$trophic.position == 1] <- "Primary Producers"
      sp.data$trophic.guild[sp.data$trophic.position >= 2 & sp.data$trophic.position < 2.3] <- "Herbivores"
      sp.data$trophic.guild[sp.data$trophic.position > 2.7 & sp.data$trophic.position < 3.3] <- "Carnivores"
      sp.data$trophic.guild[sp.data$trophic.position > 3.7 & sp.data$trophic.position < 4.3] <- "Carnivores"
      sp.data$trophic.guild[sp.data$trophic.position > 4.7 & sp.data$trophic.position < 5.3] <- "Carnivores"
      sp.data$trophic.guild[sp.data$trophic.position > 5.7 & sp.data$trophic.position < 6.3] <- "Carnivores"
      sp.data$trophic.guild[sp.data$trophic.position > 6.7 & sp.data$trophic.position < 7.3] <- "Carnivores"
      
      # table(sp.data$trophic.guild)
      
      # store results
      sim.results <- rbind(sim.results,sp.data)
      
    }# for i.replicate
    
    print(paste(date()," - niche model simulation -- richness: ", richness.levels[i.richness], ", connectance: ", connectance.levels[i.connectance], ", ",replicates," replicates completed",sep=""))
    
  }# for i.connectance
}# for i.richness

readr::write_delim(x = sim.results,path = "./results/abundances_niche_model.csv",delim = ";")
