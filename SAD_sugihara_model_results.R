#####################
# calculate hill evenness, mad, skewness for every SAD obtained with the script "SAD_cohen_abundances_model"
#####################

model.data <- readr::read_delim(file = "./results/abundances_niche_model.csv",delim = ";")

richness.levels <- unique(model.data$richness)
connectance.levels <- unique(model.data$connectance)
apportionment.levels <- unique(model.data$niche.apport)
trophic.guilds <- unique(model.data$trophic.guild)
replicates <- max(model.data$replicate)

# guild.abundances <- model.data %>% group_by(richness,connectance,trophic.guild,replicate) %>% summarise(abund = n())
# mean.abundances <- guild.abundances %>% group_by(richness,connectance,trophic.guild) %>% summarise(mean.abund = mean(abund))
# ggplot(guild.abundances) + 
#   geom_boxplot(aes(x = trophic.guild, y = abund)) + 
#   facet_grid(richness~connectance, scales="free") + 
#   NULL
#############
metrics.results <- NULL

# i.richness <- 1
# i.connectance <- 3
# i.tl <- 3
# i.rep <- 2

for(i.richness in 1:length(richness.levels)){
  for(i.connectance in 1:length(connectance.levels)){
    for(i.apport in 1:length(apportionment.levels)){
      for(i.tl in 1:length(trophic.guilds)){
        for(i.rep in 1:replicates){
          
          temp.result <- data.frame(richness.level = richness.levels[i.richness],
                                    connectance.level = connectance.levels[i.connectance],
                                    apportionment.level = apportionment.levels[i.apport],
                                    replicate = i.rep,
                                    trophic.guild = trophic.guilds[i.tl],
                                    guild.richness = 0,
                                    mad = 0,
                                    hill.evenness = 0,
                                    skewness = 0,
                                    log.mad = 0,
                                    log.hill.evenness = 0,
                                    log.skewness = 0)
          
          my.abundances <- model.data$N_predicted[model.data$richness == richness.levels[i.richness] & 
                                                    model.data$connectance == connectance.levels[i.connectance] &
                                                    model.data$niche.apport == apportionment.levels[i.apport] &
                                                    model.data$trophic.guild == trophic.guilds[i.tl] &
                                                    model.data$replicate == i.rep]
          
          # transform data for consistency
          # abundances < 0.5 -> 0
          # abundances 0.5 < x <= 1 -> 1.001
          # this way log.data can be computed and hill.diversity function does not raise errors
          my.abundances[my.abundances < 0.5] <- 0
          my.abundances[my.abundances >= 0.5 & my.abundances <= 1] <- 1.001
          
          if(sum(my.abundances)>0){
            
            temp.result$guild.richness <- sum(my.abundances)
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
              cat(i.richness,"-",i.connectance,"-",trophic.guilds[i.tl],"-",i.rep,"\n","log.data:",log.data,"\n","abundances:",my.abundances)
              stop()
            }
            temp.result$log.mad <- stats::mad(log.data,constant = 1)
            temp.result$log.skewness <- robustbase::mc(log.data)
            log.hill.diversity <- hill.diversity(log.data)
            temp.result$log.hill.evenness <- log.hill.diversity/length(my.abundances)
            
            metrics.results <- rbind(metrics.results,temp.result)
            
          }# if any abundance
        }# for i.rep
      }# for trophic.guilds[i.tl]
    }# for i.apport
  }# for i.connectance
}# for i.richness

# readr::write_delim(x = metrics.results,path = "./results/cohen_simulation_metrics_CONSTANT_DEGREE.csv",delim = ";")

readr::write_delim(x = metrics.results,path = "./results/sugihara_model_metrics.csv",delim = ";")

