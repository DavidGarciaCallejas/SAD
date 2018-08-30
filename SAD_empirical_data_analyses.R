#####################
# statistical analyses of the empirical datasets
#####################

# use log-transformed or natural abundances?
log.transformed <- FALSE

#####################

# read empirical data
site.metrics <- readr::read_delim(file = "./results/empirical_metrics.csv",
                                  delim = ";",
                                  col_types = list(col_character(),
                                                   col_character(),
                                                   col_double(),
                                                   col_double(),
                                                   col_double(),
                                                   col_double(),
                                                   col_double(),
                                                   col_double(),
                                                   col_double(),
                                                   col_double(),
                                                   col_double(),
                                                   col_double(),
                                                   col_double(),
                                                   col_double()))
site.parameters <- readr::read_delim(file = "./data/abundance datasets/all_sites.csv",delim = ";")

#####################
# combine them in a single dataset
site.metrics <- left_join(site.metrics,site.parameters)
# some sites do not have spatial or temporal extent, remove them
site.metrics <- subset(site.metrics,!is.na(dataset))
site.metrics$habitat <- as.factor(site.metrics$habitat)
site.metrics$trophic.guild <- as.factor(site.metrics$trophic.guild)
site.metrics$trophic.guild <- relevel(site.metrics$trophic.guild, ref = "plants")
site.metrics$hill.evenness.log[site.metrics$hill.evenness.log>1] <- 1

site.metrics <- droplevels(site.metrics)

# other cleaning
# 1 - guilds with more than X species, otherwise it does not make sense to compute SAD
minimum.richness <- 3
site.metrics <- subset(site.metrics,richness >= minimum.richness)

site.metrics$trophic.guild <- factor(site.metrics$trophic.guild, levels = c("plants","herbivores","omnivores","carnivores"))

if(log.transformed){
  site.metrics <- site.metrics[,c("site.ID","trophic.guild","hill.evenness.log","skewness.log","richness","spatial.extent","temporal.extent","habitat","dataset")]
  names(site.metrics)[c(3,4)] <- c("hill.evenness","skewness")
}else{
  site.metrics <- site.metrics[,c("site.ID","trophic.guild","hill.evenness","skewness","richness","spatial.extent","temporal.extent","habitat","dataset")]
}

#####################
# statistical analyses

# 0 - data exploration

# hist(site.metrics$hill.evenness, 100, col="black")
# clearly one-inflated

# test.plot1 <- ggplot(site.metrics,aes(y = hill.evenness,x = temporal.extent)) +
#   geom_point(aes(color = trophic.guild))
# test.plot1

# test.plot2 <- ggplot(no.plants.data,aes(y = hill.evenness,x = temporal.extent)) +
#   geom_point(aes(color = trophic.guild))
# test.plot2

# test.plot3 <- ggplot(no.plants.data,aes(y = hill.evenness,x = spatial.extent)) +
#   geom_point(aes(color = trophic.guild))
# test.plot3

# test.plot.4 <- ggplot(site.metrics,aes(y = hill.evenness, x= trophic.guild)) +
#   geom_boxplot(aes(group = trophic.guild))
# test.plot.4

# test.plot.41 <- ggplot(site.metrics,aes(y = hill.evenness, x= habitat)) +
#   geom_boxplot(aes(group = habitat)) + facet_grid(trophic.guild~.)
# test.plot.41
# 
# test.plot.5 <- ggplot(site.metrics,aes(y = hill.evenness, x= habitat)) +
#   geom_boxplot(aes(group = habitat))
# test.plot.5

#################
# In order to minimise the complexity of the model, it may make sense to transform the data slightly, 
# so that it lies in (0,1) and a standard beta regression can be fitted.
site.metrics$transformed.evenness <- (site.metrics$hill.evenness*(nrow(site.metrics)-1) + 0.5)/nrow(site.metrics)
# how is the transformation? does it vary much?
# hist(site.metrics$hill.evenness, 100, col="black")
# hist(site.metrics$transformed.evenness,100,col="black")
# plot(site.metrics$hill.evenness - site.metrics$transformed.evenness)
# max(site.metrics$hill.evenness - site.metrics$transformed.evenness)
# 
# # transformation looks good
# # try a beta reg model without random effects
# betareg.model <- betareg(transformed.evenness ~ trophic.guild + spatial.extent + temporal.extent + richness,data = site.metrics)
# 
# # diagnostics
# plot(betareg.model, which = 1:4, type = "pearson")
# plot(betareg.model, which = 5, type = "deviance", sub.caption = "")
# plot(betareg.model, which = 1, type = "deviance", sub.caption = "")
# summary(betareg.model)$pseudo.r.squared
# AIC(betareg.model)

# try another link function (betareg vignette suggests log-log for data with many extreme values)
# however, a test of the likelihoods says that clog-log is better
# sapply(c("logit","probit","cloglog","cauchit"),function(x) logLik(update(betareg.model, link = x)))
# sapply(c("logit", "probit", "cloglog", "cauchit", "loglog"),function(x) logLik(update(betareg.model, link = x)))

# 
# betareg.loglog <- betareg(transformed.evenness ~ trophic.guild + spatial.extent + temporal.extent + richness,data = site.metrics,link = "cloglog")
# 
# plot(betareg.loglog, which = 1:4, type = "pearson")
# plot(betareg.loglog, which = 5, type = "deviance", sub.caption = "")
# plot(betareg.loglog, which = 1, type = "deviance", sub.caption = "")
# summary(betareg.loglog)$pseudo.r.squared
# AIC(betareg.loglog)

# try a more detailed formula for the precision parameter
# betareg.loglog.pre <- betareg(transformed.evenness ~ trophic.guild + spatial.extent + temporal.extent + richness | trophic.guild + spatial.extent + temporal.extent + richness,
#                               data = site.metrics,link = "cloglog")
# 
# plot(betareg.loglog.pre, which = 1:4, type = "pearson")
# plot(betareg.loglog.pre, which = 5, type = "deviance", sub.caption = "")
# plot(betareg.loglog.pre, which = 1, type = "deviance", sub.caption = "")
# summary(betareg.loglog.pre)$pseudo.r.squared
# AIC(betareg.loglog.pre)
# evenness.tl.emm.betareg <- emmeans(betareg.loglog.pre, "trophic.guild")
# pairs(evenness.tl.emm.betareg)

# also, try another link functions for the precision parameter
# sapply(c("identity","log","sqrt"),function(x) logLik(update(betareg.model, link.phi = x)))

# best one seems the last one
# summary(betareg.loglog.pre)

####################
# now, define the full model, see main text for details about the mu and sigma parameters.

evenness.beta <- gamlss(formula = transformed.evenness ~ trophic.guild + richness + spatial.extent + temporal.extent, 
                        sigma.formula = transformed.evenness ~ trophic.guild + spatial.extent + temporal.extent + richness,
                        family = BE(mu.link = cloglog),data = site.metrics)

# with interaction trophic.guild * richness??
# evenness.beta.interaction <- gamlss(formula = transformed.evenness ~ trophic.guild*richness + spatial.extent + temporal.extent, #+ random(habitat),
#                         sigma.formula = transformed.evenness ~ trophic.guild + spatial.extent + temporal.extent + richness,# + random(habitat),
#                         family = BE(mu.link = cloglog),data = site.metrics)
# 
# AIC(evenness.beta.interaction,evenness.beta)

# including the interaction only marginally improves the model,
# so I use the simpler model
summary(evenness.beta)
Rsq(evenness.beta)

evenness.tl.emm <- emmeans(evenness.beta, "trophic.guild")
pairs(evenness.tl.emm)

#################
# model comparison by AIC with some null models
evenness.null.1 <- gamlss(formula = transformed.evenness ~ 1, 
                          family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.2 <- gamlss(formula = transformed.evenness ~ trophic.guild, 
                        sigma.formula = transformed.evenness ~ trophic.guild,
                        family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.3 <- gamlss(formula = transformed.evenness ~ richness, 
                          sigma.formula = transformed.evenness ~ richness,
                          family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.4 <- gamlss(formula = transformed.evenness ~ spatial.extent, 
                          sigma.formula = transformed.evenness ~ temporal.extent,
                          family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.5 <- gamlss(formula = transformed.evenness ~ temporal.extent, 
                          sigma.formula = transformed.evenness ~ temporal.extent,
                          family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.6 <- gamlss(formula = transformed.evenness ~ trophic.guild + richness, 
                          sigma.formula = transformed.evenness ~ trophic.guild + richness,
                          family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.7 <- gamlss(formula = transformed.evenness ~ trophic.guild + spatial.extent, 
                          sigma.formula = transformed.evenness ~ trophic.guild + spatial.extent,
                          family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.8 <- gamlss(formula = transformed.evenness ~ trophic.guild + temporal.extent, 
                          sigma.formula = transformed.evenness ~ trophic.guild + temporal.extent,
                          family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.9 <- gamlss(formula = transformed.evenness ~ richness + spatial.extent, 
                          sigma.formula = transformed.evenness ~ richness + spatial.extent,
                          family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.10 <- gamlss(formula = transformed.evenness ~ richness + temporal.extent, 
                          sigma.formula = transformed.evenness ~ richness + temporal.extent,
                          family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.11 <- gamlss(formula = transformed.evenness ~ spatial.extent + temporal.extent, 
                          sigma.formula = transformed.evenness ~ spatial.extent + temporal.extent,
                          family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.12 <- gamlss(formula = transformed.evenness ~ trophic.guild + richness + spatial.extent, 
                          sigma.formula = transformed.evenness ~ trophic.guild + richness + spatial.extent,
                          family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.13 <- gamlss(formula = transformed.evenness ~ trophic.guild + richness + temporal.extent,
                          sigma.formula = transformed.evenness ~ trophic.guild + richness + temporal.extent,
                          family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.14 <- gamlss(formula = transformed.evenness ~ trophic.guild + spatial.extent + temporal.extent, 
                          sigma.formula = transformed.evenness ~ trophic.guild + spatial.extent + temporal.extent,
                          family = BE(mu.link = cloglog),data = site.metrics)
evenness.null.15 <- gamlss(formula = transformed.evenness ~ richness + spatial.extent + temporal.extent, 
                          sigma.formula = transformed.evenness ~ richness + spatial.extent + temporal.extent,
                          family = BE(mu.link = cloglog),data = site.metrics)

AIC(evenness.beta,evenness.null.1,evenness.null.2,evenness.null.3,evenness.null.4,evenness.null.5,
    evenness.null.6,evenness.null.7,evenness.null.8,evenness.null.9,evenness.null.10,evenness.null.11,
    evenness.null.12,evenness.null.13,evenness.null.14,evenness.null.15)

#################
#################
# skewness
hist(site.metrics$skewness,100,col="black")
# hist(LinMap(site.metrics$skewness,0,1),100)

# 
# test.plot1 <- ggplot(site.metrics,aes(y = skewness.log,x = temporal.extent)) +
#   geom_point(aes(color = trophic.guild))
# test.plot1
# 
# test.plot3 <- ggplot(site.metrics,aes(y = skewness.log,x = spatial.extent)) +
#   geom_point(aes(color = trophic.guild))
# test.plot3
# 
# test.plot.4 <- ggplot(site.metrics,aes(y = skewness.log, x= trophic.guild)) +
#   geom_boxplot(aes(group = trophic.guild))
# test.plot.4
# 
# test.plot.41 <- ggplot(site.metrics,aes(y = skewness.log, x= trophic.guild)) +
#   geom_boxplot(aes(group = trophic.guild)) + facet_wrap(~habitat)
# test.plot.41
# 
# test.plot.5 <- ggplot(site.metrics,aes(y = skewness.log, x= habitat)) +
#   geom_boxplot(aes(group = habitat))
# test.plot.5
#################
# clearly three modes. Divide in three categories
site.metrics$skewness.categories <- character(nrow(site.metrics))

for(i.row in 1:nrow(site.metrics)){
  if(site.metrics$skewness[i.row] < -0.5){
    site.metrics$skewness.categories[i.row] <- "[-1,-0.5)"
  }else if(site.metrics$skewness[i.row] >= -0.5 & site.metrics$skewness[i.row] <= 0.5){
    site.metrics$skewness.categories[i.row] <- "[-0.5,0.5]"
  }else if(site.metrics$skewness[i.row] > 0.5){
    site.metrics$skewness.categories[i.row] <- "(0.5,1]"
  }
}# for

#################
# multinomial logistic regression with the same predictors as evenness model
site.metrics$skewness.categories <- as.factor(site.metrics$skewness.categories)
site.metrics$skewness.categories <- relevel(site.metrics$skewness.categories,ref = "[-0.5,0.5]")

sk.multinom <- multinom(skewness.categories ~ trophic.guild + spatial.extent + temporal.extent + richness, 
                        data = site.metrics,model=TRUE)
summary(sk.multinom)

##########
# plot these categories
sk.cat.data <- site.metrics %>% group_by(trophic.guild,skewness.categories) %>% summarise(num.sites = n())
sk.cat.data$skewness.categories <- factor(sk.cat.data$skewness.categories, levels = c("[-1,-0.5)","[-0.5,0.5]","(0.5,1]"))
my.palette <- c("darkgreen","#009E73","#E69F00","#D55E00")

sk.cat.plot <- ggplot(sk.cat.data) + 
  
  geom_col(aes(x = trophic.guild, y = num.sites, fill = trophic.guild)) +
  facet_wrap(~skewness.categories) +
  
  # geom_col(aes(x = skewness.categories, y = num.sites, fill = trophic.guild)) +
  # facet_wrap(~trophic.guld, scales = "free_y") +
  
  scale_fill_manual(values = my.palette)+
  xlab("skewness levels") + ylab("number of sites") +
  #DGC::theme_Publication() +
  scale_x_discrete(breaks=NULL)+
  scale_y_continuous(expand = c(0.01, 0)) +
  theme(strip.background = element_blank())+
  NULL
sk.cat.plot

tiff("./results/images/skewness_categories.tiff", res=600, compression = "lzw", width = 4500, height = 3500, units = "px")
sk.cat.plot
dev.off()

#########################################
# AIC comparison against some null models

skewness.null.1 <- multinom(skewness.categories ~ 1, data = site.metrics,model=TRUE)#,Hess = TRUE)
skewness.null.2 <- multinom(skewness.categories ~ trophic.guild,data = site.metrics,model = TRUE)
skewness.null.3 <- multinom(skewness.categories ~ richness, data = site.metrics,model = TRUE)
skewness.null.4 <- multinom(skewness.categories ~ spatial.extent, data = site.metrics,model = TRUE)
skewness.null.5 <- multinom(skewness.categories ~ temporal.extent, data = site.metrics,model = TRUE)
skewness.null.6 <- multinom(skewness.categories ~ trophic.guild + richness, data = site.metrics,model = TRUE)
skewness.null.7 <- multinom(skewness.categories ~ trophic.guild + spatial.extent, data = site.metrics,model = TRUE)
skewness.null.8 <- multinom(skewness.categories ~ trophic.guild + temporal.extent, data = site.metrics,model = TRUE)
skewness.null.9 <- multinom(skewness.categories ~ richness + spatial.extent, data = site.metrics,model = TRUE)
skewness.null.10 <- multinom(skewness.categories ~ richness + temporal.extent, data = site.metrics,model = TRUE)
skewness.null.11 <- multinom(skewness.categories ~ spatial.extent + temporal.extent, data = site.metrics,model = TRUE)
skewness.null.12 <- multinom(skewness.categories ~ trophic.guild + richness + spatial.extent, data = site.metrics,model = TRUE)
skewness.null.13 <- multinom(skewness.categories ~ trophic.guild + temporal.extent + richness, data = site.metrics,model = TRUE)
skewness.null.14 <- multinom(skewness.categories ~ trophic.guild + spatial.extent + temporal.extent, data = site.metrics,model = TRUE)
skewness.null.15 <- multinom(skewness.categories ~ richness + spatial.extent + temporal.extent, data = site.metrics,model = TRUE)


CatgoF::goFit(sk.multinom,skewness.null.1,skewness.null.2,skewness.null.3,skewness.null.4,skewness.null.5,skewness.null.6,skewness.null.7,skewness.null.8,skewness.null.9,skewness.null.10,skewness.null.11,
              skewness.null.12,skewness.null.13,skewness.null.14,skewness.null.15,
              crit = "AIC",display = "all")

# best model according to AIC is sk.null.6, trophic.guild + richness 
summary(skewness.null.6)

z <- summary(skewness.null.6)$coefficients/summary(skewness.null.6)$standard.errors
# 2-tailed Wald z tests to test significance of coefficients
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

nnet.mod.loglik <- nnet:::logLik.multinom(skewness.null.6)
nnet.mod0 <- multinom(skewness.categories ~ 1, site.metrics)
nnet.mod0.loglik <- nnet:::logLik.multinom(nnet.mod0)
r_squared_multinom <- as.numeric(1 - nnet.mod.loglik/nnet.mod0.loglik)

# report results from summary, alongside z-score and p-values

################
# trophic guild effects

# this is the appropriate comparison
sk.comp <- emmeans(skewness.null.6, ~  trophic.guild | skewness.categories, mode = "prob")
pairs(sk.comp)
