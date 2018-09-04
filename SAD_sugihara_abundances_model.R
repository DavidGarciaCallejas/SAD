
# theoretical distributions of abundances in increasing trophic levels
# given body sizes of all levels, and abundances of the basal level

#####################
#####################
# auxiliary functions

#####################
# the niche model for food web structure


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

connectance <- 0.2
S <- 10
shape.weibull.min <- 0.15
shape.weibull.max <- 0.2
scale.weibull <- 4.7 * 100

tt <- NicheModel(S,connectance)
interaction.matrix <- tt[[1]]
niche <- tt[[2]]

abundances <- data.frame(ID = 1:S,M = niche,abundance = 0)

basal <- which(colSums(interaction.matrix) == 0)
abundances$abundance[basal] <- rweibull(length(basal),
                                        shape = runif(1,shape.weibull.min,shape.weibull.max),
                                        scale = scale.weibull)

for(i.basal in 1:length(basal)){
  # predator.set <- which(interaction.matrix[i.basal,] == 1)
  predator.set <- data.frame(predator.id = which(interaction.matrix[i.basal,] == 1), specialization = 0)
  if(nrow(predator.set)>1){
    # resource is partitioned according to the degree of specialization
    # and this is given by the number of links of the predators
    for(i.predator in 1:nrow(predator.set)){
      predator.set$specialization[i.predator] <- sum(interaction.matrix[,predator.set$predator.id[i.predator]] == 1)
    }# for each predator
    
    resource.available <- abundances$M[i.basal] * abundances$abundance[i.basal]
    # dominance preemption, sugihara 1980, tokeshi 1990
    resource.partition <- nicheApport::dominancePreemp(resource.available,nrow(predator.set))
    # if i ever use dominance decay, beware that it fails for n=2
    predator.set <- arrange(predator.set,specialization)
    predator.set$assigned <- resource.partition
    
    abundances$abundance[abundances$ID %in% predator.set$predator.id] 
    
    # # those with assigned resources
    # predator.set$assigned <- FALSE
    # # how much of the resource still available
    # available.resource <- 1
    # # there are (predator-1) breaks
    # for(i.break in 1:(length(predator.set)-1)){
    #   
    #   # dominant species will take the biggest slice
    #   dominant <- predator.set$predator.id[predator.set$specialization == min(predator.set$specialization) & predator.set$assigned == FALSE]
    #   
    #   # if there are two or more equally dominant, randomize
    #   if(length(dominant)>1){
    #     dominant <- sample(dominant,1)
    #   }
    #   
    #   resource <- abundances$M[i.basal] * abundances$abundance[i.basal] * resource.available
    #   
    #   dominant.proportion <- runif(1,0.5,0.99)
    #   
    #   # which break is to be split after this break? 
    #   # 1 is the dominant
    #   my.break <- sample(1:2,1)
    #   
    #   # if the dominant piece is to be split after this break, 
    #   if(my.break == 1){
    #     
    #   }else{
    #     abundances$abundance[abundances$ID == dominant] <- abundances$abundance[abundances$ID == dominant] + dominant.proportion*resource
    #   }
    #   
    #   
    #   
    # }
      
  }else if(nrow(predator.set) == 1){
    abundances$abundance[abundances$ID == predator.set] <- abundances$abundance[i.basal]
  }else if(length(predator.set) == 0){
    # no predators, so nothing to do
  }
  
}# for i.basal

#####################
# this function generates samples from a truncated poisson distribution
# N = number of samples
# lambda = poisson rate parameter
# k = lower limit (right-truncated)

rtpois <- function(N, lambda, k) qpois(runif(N, ppois(k, lambda), 1), lambda)

#####################
#####################
# simulations for obtaining SADs based on Cohen's formula
# communities with:
# 4 trophic levels
# varying degrees of (1) vulnerability and (2) body mass-trophic level relationships

# qualitatively similar to fig. 1 of turney and buddle 2016
richness <- 100
richness.tl.proportion <- c(.35,.25,.25,.15)

# body size levels
# see Riede et al. (2011)
# 1 - body mass increases with trophic level
# 2 - predator-prey mass ratio generally decreases with trophic level
# such that body size of predators and their prey are more similar up in the food chain

# if p = predator/prey mass ratio,
# log(p,10) = a + b*trophic.level
# b is the slope of the relationship.
# calculate b from table 2 of Riede et al (2011)
mass.ratio.exponent <- mean(c(-0.82,-0.89,-0.87,-0.49,-0.62,-0.71,-1.5,-0.39,-0.49,0.56))
# the exact form will be calculated at each timestep, given stochastic values of the intercept
# body size of basal trophic level
basal.mass <- 1

# how many replicates of each configuration
replicates <- 100

# vulnerability levels
# min: one predator per prey
# max: all predators prey upon all preys
# in between, there are three levels, with various percentages of vulnerability
# these three percentages are given by "vulnerability.connectances"
vulnerability.levels <- c("min","low","intermediate","high","max")
vulnerability.connectances <- c(0.1,0.2,0.5)

# vulnerability, constant for the whole trophic level or varying degree, depending on abundance?
constant.vulnerability <- FALSE

# in this version, only increasing body size, following Riede et al. (2011)
bodysize.tl.levels <- c("increasing")

# which distribution for generating basal abundances?
basal.dist <- "weibull"

# these parameters are taken from "SAD_cohen_fit_basal_dist.R"
# with the additional point that I multiplied the scale parameter by an arbitrary value
# in order to get maximum basal abundances on the range of 10^3
shape.weibull.min <- 0.15
shape.weibull.max <- 0.2
scale.weibull <- 4.7 * 100

# if using exponential distribution...
# rate of the exponential distribution for generating basal abundances
exp.rate <- 0.005

#####################
#####################

sim.results <- NULL

i.vul <- 2
i.size <- 1
i.rep <- 1

for(i.vul in 1:length(vulnerability.levels)){
  for(i.size in 1:length(bodysize.tl.levels)){
    
    for(i.rep in 1:replicates){
      
      # this is arbitrary
      mass.ratio.intercept <- runif(1,5,10)
      
      # calculate body mass of higher trophic levels
      size.cat <- c(basal.mass,0,0,0)
      for(i.tl in 2:4){
        my.ratio <- mass.ratio.intercept * i.tl^mass.ratio.exponent
        my.cat <- size.cat[i.tl-1]*my.ratio
        size.cat[i.tl] <- my.cat
      }
      
      
      # species info
      sp.data <- data.frame(vulnerability.level = vulnerability.levels[i.vul],
                            body.size.tl = bodysize.tl.levels[i.size],
                            replicate = i.rep,
                            ID = 1:richness,trophic.level = c(rep(1,richness * richness.tl.proportion[1]),
                                                              rep(2,richness * richness.tl.proportion[2]),
                                                              rep(3,richness * richness.tl.proportion[3]),
                                                              rep(4,richness * richness.tl.proportion[4])),M = 0, V = 0, N_predicted = 0)
      # body mass - trophic level relationship
      # allowing some variability (0.1*mean)
      sp.data$M[sp.data$trophic.level == 1] <- rnorm(nrow(sp.data[sp.data$trophic.level == 1,]),size.cat[1],size.cat[1]*.1)
      sp.data$M[sp.data$trophic.level == 2] <- rnorm(nrow(sp.data[sp.data$trophic.level == 2,]),size.cat[2],size.cat[2]*.1)
      sp.data$M[sp.data$trophic.level == 3] <- rnorm(nrow(sp.data[sp.data$trophic.level == 3,]),size.cat[3],size.cat[3]*.1)
      sp.data$M[sp.data$trophic.level == 4] <- rnorm(nrow(sp.data[sp.data$trophic.level == 4,]),size.cat[4],size.cat[4]*.1)    
      
      # abundance of the basal trophic level
      if(basal.dist == "weibull"){
        sp.data$N_predicted[sp.data$trophic.level == 1] <- rweibull(nrow(sp.data[sp.data$trophic.level == 1,]),
                                                                    shape = runif(1,shape.weibull.min,shape.weibull.max),
                                                                    scale = scale.weibull)
      }else if(basal.dist == "exp"){
        sp.data$N_predicted[sp.data$trophic.level == 1] <- rexp(nrow(sp.data[sp.data$trophic.level == 1,]),exp.rate)
      }
      
      # generate links according to the vulnerability level
      interaction.matrix <- matrix(0,nrow = richness,ncol = richness)
      
      # basal level does not have prey
      for(i.tl in 2:4){
        prey.tl <- subset(sp.data,trophic.level == i.tl-1)
        
        # species in the upper tl/ or, in general, in a higher tl if omnivory!!
        # predator.set <- sp.data$ID[sp.data$trophic.level == (i.tl)]
        predator.set <- sp.data$ID[sp.data$trophic.level > i.tl]
        
        # average number of links
        if(vulnerability.levels[i.vul] == "min"){
          my.size <- 1
        }else if(vulnerability.levels[i.vul] == "low"){
          my.size <- round(length(predator.set) * vulnerability.connectances[1])
        }else if(vulnerability.levels[i.vul] == "intermediate"){
          my.size <- round(length(predator.set) * vulnerability.connectances[2])
        }else if(vulnerability.levels[i.vul] == "high"){
          my.size <- round(length(predator.set) * vulnerability.connectances[3])
        }else{
          my.size <- length(predator.set)
        }
        # in case it is rounded to 0
        my.size <- ifelse(my.size==0,1,my.size)
        
        # same vulnerability (degree) for all species?
        # if not, generate random samples from a truncated poisson distribution and sort
        # by abundances
        if(constant.vulnerability){
          prey.tl$V <- rep(my.size,nrow(prey.tl))
        }else{
          prey.tl.links <- rtpois(N = nrow(prey.tl),lambda = my.size,k = 0)
          prey.tl.links[prey.tl.links > length(predator.set)] <- length(predator.set)
          prey.tl <- arrange(prey.tl,desc(N_predicted))
          prey.tl$V <- sort(prey.tl.links,decreasing = TRUE)
          
          prey.tl <- arrange(prey.tl,ID)
        }
        
        # now, fill the interactions of each species
        for(i.sp in 1:nrow(prey.tl)){
          
          if(length(predator.set)>1){
            unconnected.predators <- predator.set[which(colSums(interaction.matrix[,predator.set]) >= 0)]
          }else{
            unconnected.predators <- ifelse(sum(interaction.matrix[,predator.set])<1,predator.set,0)
          }
          
          # if number of links is less than the number of unconnected predators, link them
          # otherwise, choose from the overall predator pool.
          # don't bother too much.
          if(prey.tl$V[i.sp] <= length(unconnected.predators)){
            my.predator <- sample.vec(x = unconnected.predators,size = prey.tl$V[i.sp],replace=F)
          }else{
            my.predator <- sample.vec(x = predator.set,size = prey.tl$V[i.sp],replace=F)
          }
          
          interaction.matrix[i.sp,my.predator] <- -1
          interaction.matrix[my.predator,i.sp] <- 1
          
          sp.data$V[sp.data$ID == prey.tl$ID[i.sp]] <- prey.tl$V[i.sp]
          
        }# for each species
        
        # calculate niche apportionment and abundances of this trophic level
        # sugihara 1980 is the main reference.
        
        for(i.sp in 1:nrow(sp.data)){
          if(sp.data$trophic.level[i.sp] == i.tl){
            # compute expected abundances...
            # prey species
            prey.sp <- which(interaction.matrix[,i.sp] == -1)
            if(length(prey.sp) != 0){
              prey.data <- sp.data[sp.data$ID %in% prey.sp,]
              sp.data$N_predicted[i.sp] <- CohenAbundance(species.mass = sp.data$M[i.sp],prey.data = prey.data,type = "N")
            }# if-else prey
            
          }# if my trophic level
        }# for i.sp
        
      }# for each trophic level
      
      # store results
      sim.results <- rbind(sim.results,sp.data)
      
    }# for i.rep
  }# for i.size
}# for i.vul

if(constant.vulnerability){
  readr::write_delim(x = sim.results,path = "./results/cohen_simulation_abundances_CONSTANT_DEGREE.csv",delim = ";")
}else{
  readr::write_delim(x = sim.results,path = "./results/cohen_simulation_abundances_VARIABLE_DEGREE.csv",delim = ";")
}






