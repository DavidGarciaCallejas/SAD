###################
# evaluate goodness of fit of cohen abundances with the tuesday lake data
###################

# alpha = 1 or alpha = 0.75?
alpha.metabolism <- 1

if(alpha.metabolism == 1){
  my.file.name <- "./results/cohen_tuesday_lake_abundances.csv"
}else{
  my.file.name <- "./results/cohen_tuesday_lake_abundances_with_alpha.csv"
}

predicted.abundances <- readr::read_delim(file = my.file.name,delim = ";")

##############################
##############################

# hist(log(predicted.abundances$N_pred_full),100,col="black")

fixed.model <- lm(log(N_pred_full) ~ log(N_obs) + dataset,data = predicted.abundances)

# summary(fixed.model)

#################
# plot the fitted model
predicted.abundances$fit <- predict(fixed.model)

residuals.lm <- ggplot(predicted.abundances,aes(x = log(N_obs), y = log(N_pred_full), color = dataset))

residuals.lm <- residuals.lm + 
  geom_point() +
  geom_line(aes(y = fit)) + 
  # geom_smooth(method = "lm",aes(color = dataset),alpha = 0.3) +
  geom_abline(intercept=0, slope=1, color = "grey70",linetype = "dashed") +
  # facet_wrap(~dataset) + 
  scale_color_manual(name="",labels = c("1984","1986"),values = c("#0072B2","#D55E00")) + 
  xlab("log(observed abundances)") + ylab("log(predicted abundances)") + 
  #guides(color = F) + 
  #DGC::theme_Publication()+
  NULL

my.plot.name <- ifelse(alpha.metabolism == 1,
                       "./results/images/Cohen_fit_lm.tiff",
                       "./results/images/Cohen_fit_lm_with_alpha.tiff")

tiff(my.plot.name, res=600, compression = "lzw", width = 3000, height = 2500, units = "px")
residuals.lm
dev.off()
