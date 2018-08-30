#####################
# calculate hill evenness, mad, skewness for every SAD obtained with the script "SAD_cohen_abundances_model"
#####################

model.data <- readr::read_delim(file = "./results/cohen_simulation_abundances_VARIABLE_DEGREE.csv",delim = ";")

vulnerability.levels <- unique(model.data$vulnerability.level)
bodysize.tl.levels <- unique(model.data$body.size.tl)
replicates <- max(model.data$replicate)
trophic.levels <- max(model.data$trophic.level)

#############
metrics.results <- NULL

# i.vul <- 1
# i.size <- 3
# i.tl <- 3
# i.rep <- 2

for(i.vul in 1:length(vulnerability.levels)){
  for(i.size in 1:length(bodysize.tl.levels)){
    for(i.tl in 1:trophic.levels){
      for(i.rep in 1:replicates){
        
        temp.result <- data.frame(vulnerability.level = vulnerability.levels[i.vul],
                                  body.size.tl = bodysize.tl.levels[i.size],
                                  trophic.level = i.tl,
                                  replicate = i.rep,
                                  mad = 0,
                                  hill.evenness = 0,
                                  skewness = 0,
                                  log.mad = 0,
                                  log.hill.evenness = 0,
                                  log.skewness = 0)
        
        my.abundances <- model.data$N_predicted[model.data$vulnerability.level == vulnerability.levels[i.vul] & 
                                                  model.data$body.size.tl == bodysize.tl.levels[i.size] & 
                                                  model.data$trophic.level == i.tl &
                                                  model.data$replicate == i.rep]
        
        # transform data for consistency
        # abundances < 0.5 -> 0
        # abundances 0.5 < x <= 1 -> 1.001
        # this way log.data can be computed and hill.diversity function does not raise errors
        my.abundances[my.abundances < 0.5] <- 0
        my.abundances[my.abundances >= 0.5 & my.abundances <= 1] <- 1.001
        
        if(sum(my.abundances)>0){
        
        temp.result$mad <- stats::mad(my.abundances,constant = 1)
        temp.result$skewness <- robustbase::mc(my.abundances)
        my.hill.diversity <- hill.diversity(my.abundances)
        temp.result$hill.evenness <- my.hill.diversity/length(my.abundances)
        
        # metrics for log data:
        log.data <- log(my.abundances[my.abundances != 0],2)
        # in case abundance exactly 1, add a small increment so it is not taken as 0
        # it shouldn't happen anyway with the above transformation
        log.data[log.data == 0] <- 0.001  
        if(sum(is.na(log.data))>0){
          cat(i.vul,"-",i.size,"-",i.tl,"-",i.rep,"\n","log.data:",log.data,"\n","abundances:",my.abundances)
          stop()
        }
        temp.result$log.mad <- stats::mad(log.data,constant = 1)
        temp.result$log.skewness <- robustbase::mc(log.data)
        log.hill.diversity <- hill.diversity(log.data)
        temp.result$log.hill.evenness <- log.hill.diversity/length(my.abundances)
        
        metrics.results <- rbind(metrics.results,temp.result)
        
        }# if any abundance
      }# for i.rep
    }# for i.tl
  }# for i.size
}# for i.vul

# readr::write_delim(x = metrics.results,path = "./results/cohen_simulation_metrics_CONSTANT_DEGREE.csv",delim = ";")

readr::write_delim(x = metrics.results,path = "./results/cohen_simulation_metrics_VARIABLE_DEGREE.csv",delim = ";")

