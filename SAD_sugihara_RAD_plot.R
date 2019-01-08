#####################
# calculate rank-abundance plots for SADs obtained with the script "SAD_sugihara_abundances_model"
#####################

d1 <- readr::read_delim(file = "./results/abundances_niche_model_trophic_guild_DD.csv",delim = ";")
d2 <- readr::read_delim(file = "./results/abundances_niche_model_trophic_guild_DP.csv",delim = ";")
d3 <- readr::read_delim(file = "./results/abundances_niche_model_trophic_guild_RF.csv",delim = ";")
model.data <- bind_rows(d1,d2,d3)

richness.levels <- unique(model.data$richness)
connectance.levels <- unique(model.data$connectance)
apportionment.levels <- unique(model.data$niche.apport)
trophic.guilds <- unique(model.data$trophic.guild)
replicates <- max(model.data$replicate)

rank.abundances <- model.data[,c("richness","connectance","niche.apport","replicate","ID","N_predicted","trophic.guild")]
names(rank.abundances)[6] <- "abundance"

# for averaging
average.abundances <- NULL

for(i.richness in 1:length(richness.levels)){
  for(i.connectance in 1:length(connectance.levels)){
    for(i.apport in 1:length(apportionment.levels)){
      for(i.tl in 1:length(trophic.guilds)){
        
        my.data <- subset(rank.abundances,richness == richness.levels[i.richness] & 
                            connectance == connectance.levels[i.connectance] &
                            niche.apport == apportionment.levels[i.apport] &
                            trophic.guild == trophic.guilds[i.tl])
        
        max.sp <- richness.levels[i.richness]
        avg.abund <- data.frame(abundance = numeric(max.sp))
        
        for(i.rep in 1:max(my.data$replicate)){
          temp.abund <- my.data$abundance[my.data$replicate == i.rep]
          temp.abund <- sort(temp.abund,decreasing = T)
          if(length(temp.abund) < max.sp){
            temp.abund[(length(temp.abund)):max.sp] <- 0
          }
          avg.abund$abundance <- avg.abund$abundance + temp.abund
        }
        avg.abund$richness <- richness.levels[i.richness]
        avg.abund$connectance <- connectance.levels[i.connectance]
        avg.abund$niche.apport <- apportionment.levels[i.apport]
        avg.abund$trophic.guild <- trophic.guilds[i.tl]
        
        avg.abund$abundance <- avg.abund$abundance/max(my.data$replicate)
        avg.abund <- subset(avg.abund,abundance > 0)
        
        average.abundances <- bind_rows(average.abundances,avg.abund)
      }
    }
  }
}

average.abundances <- average.abundances[,c("richness","connectance","niche.apport","trophic.guild","abundance")]


rank.abundances <- average.abundances %>% group_by(richness,connectance,niche.apport,trophic.guild) %>% mutate(relative.abund = abundance/sum(abundance))
rank.abundances <- rank.abundances %>% group_by(richness,connectance,niche.apport,trophic.guild) %>% mutate(species.rank = rank(-relative.abund,ties.method = "first"))
rank.abundances <- arrange(rank.abundances,connectance,niche.apport,trophic.guild,desc(relative.abund))

selected.rad.plot <- ggplot(rank.abundances,aes(x = species.rank,y = relative.abund, group = niche.apport)) +
  geom_line(aes(color = niche.apport), size = 1.1) +#(aes(linetype = site.ID)) +
  geom_point(aes(fill = niche.apport),shape = 21, size = 1.5) +#(aes(shape = site.ID)) + 
  facet_grid(trophic.guild~connectance+richness,scales = "free")+
  scale_color_manual(values = my.palette)+#, labels = c("plants", "herbivores", "omnivores", "carnivores (inv)", "carnivores (vert)"))+
  scale_fill_manual(values = my.palette)+#, labels = c("plants", "foli/granivores", "frugi/nectarivores", "omnivores", "carnivores (inv)", "carnivores (vert)"))+
  #c("primary\nproducers", "plant/seed\neaters", "frugivores", "omnivores", "carnivores\n(invertebrates)", "carnivores\n(vertebrates)"))+
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
  xlab("species rank") + ylab("relative abundance") +
  #DGC::theme_Publication()+
  theme(strip.background = element_blank())+#,strip.text.x = element_blank()) +
  #guides(color=FALSE)+#, fill = FALSE)+
  NULL

selected.rad.plot

# prepare data

# rank.abundances <- model.data[,c("richness","connectance","niche.apport","replicate","ID","N_predicted","trophic.guild")]
# 
# rank.abundances <- rank.abundances %>% group_by(richness,connectance,niche.apport,replicate,trophic.guild) %>% mutate(relative.abund = abundance/sum(abundance))
# rank.abundances <- rank.abundances %>% group_by(site.ID,trophic.guild) %>% mutate(species.rank = rank(-relative.abund,ties.method = "first"))

#############
metrics.results <- NULL

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

readr::write_delim(x = metrics.results,path = "./results/sugihara_model_metrics_DD.csv",delim = ";")
# readr::write_delim(x = metrics.results,path = "./results/sugihara_model_metrics_DP.csv",delim = ";")
# readr::write_delim(x = metrics.results,path = "./results/sugihara_model_metrics_RF.csv",delim = ";")

