
metrics.results <- readr::read_delim(file = "./results/sugihara_model_metrics.csv",delim = ";")

richness.levels <- unique(metrics.results$richness.level)
connectance.levels <- unique(metrics.results$connectance.level)
apportionment.levels <- unique(metrics.results$apportionment.level)
replicates <- max(metrics.results$replicate)
trophic.guilds <- unique(metrics.results$trophic.guild)

metrics.results$richness.level <- factor(metrics.results$richness.level)
metrics.results$connectance.level <- factor(metrics.results$connectance.level)
metrics.results$apportionment.level <- factor(metrics.results$apportionment.level)
metrics.results$trophic.guild <- factor(metrics.results$trophic.guild, levels = c("primary Producers","herbivores","omnivores","carnivores"))
metrics.results$trophic.guild <- plyr::revalue(metrics.results$trophic.guild,c("primary Producers" = "primary producers"))

#####

metrics.model.data <- droplevels(subset(metrics.results,trophic.guild != "primary producers"))

# null.model <- lmer(hill.evenness ~ vulnerability.level + (1 | replicate), data = metrics.model.data, REML = FALSE)
# tl.model <- lmer(hill.evenness ~ vulnerability.level + trophic.level + (1 | replicate), data = metrics.model.data, REML = FALSE)
# 
# anova(null.model,tl.model)

evenness.sugihara <- lmerTest::lmer(hill.evenness ~ richness.level + connectance.level + trophic.guild + (1 | replicate), data = metrics.model.data, REML = FALSE)
skewness.sugihara <- lmerTest::lmer(skewness ~  richness.level + connectance.level + trophic.guild + (1 | replicate), data = metrics.model.data, REML = FALSE)

# MuMIn::r.squaredGLMM(evenness.sugihara)
# MuMIn::r.squaredGLMM(skewness.sugihara)


