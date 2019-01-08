
# S <- 10
# connectance <- 0.5
# niche <- sort(runif(S,0,1))
# min.primary <- .2

NicheModel <- function(S,connectance,min.primary = 1,niche = NULL){
  
  interaction.matrix <- matrix(0,nrow=S,ncol=S)
  
  # niche axis
  if(is.null(niche)){
    my.niche <- sort(runif(S,0,1))
  }else{
    my.niche <- sort(niche)
  }  
  
  # Range
  beta = 1/(2*connectance)-1
  range = rbeta(S,1,beta)*my.niche
  
  #Centroid
  centroid = runif(S,range/2,my.niche)
  
  # Make the first niche a producer
  if(is.integer(min.primary)){
    range[1:min.primary] = 0.000000001
  }else {
    range[1:round(min.primary*S)] = .0000001
  }
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

S <- 100
connectance <- 0.5
niche <- sort(runif(S,0,1))
min.primary <- .2

interaction.matrix <- NicheModel(S = S,connectance = connectance,min.primary = min.primary,niche = niche)
primary.prod <- sum(colSums(interaction.matrix) == 0)
potential.links = (S^2 - S)
effective.connectance = sum(interaction.matrix)/potential.links



