
# theoretical distributions of abundances in increasing trophic levels
# given body sizes of all levels, and abundances of the basal level

#####################
#####################
# auxiliary functions

# this function generates samples from a truncated poisson distribution
# N = number of samples
# lambda = poisson rate parameter
# k = lower limit (right-truncated)

rtpois <- function(N, lambda, k) qpois(runif(N, ppois(k, lambda), 1), lambda)

#####################
#####################

CohenAbundance <- function(species.mass,prey.data,type = "N"){
  
  if(nrow(prey.data)>0){
    sum.biomass <- sum((prey.data$N*prey.data$M)/prey.data$V)
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
        predator.set <- sp.data$ID[sp.data$trophic.level == (i.tl)]
        # predator.set <- sp.data$ID[sp.data$trophic.level > prey.tl]
        
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
        
        # now, fill each species
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
        
        # calculate cohen abundances of this trophic level
        
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






