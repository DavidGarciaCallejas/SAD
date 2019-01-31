skewed.resources <- readr::read_delim(file = "./results/abundances_niche_model_DP_DD_RF_RA.csv",delim = ";")
uniform.resources <- readr::read_delim(file = "./results/abundances_niche_model_DP_DD_RF_RA_uniform.csv",delim = ";")

skewed.resources$resource.distribution <- "skewed"
uniform.resources$resource.distribution <- "uniform"

all.data <- bind_rows(uniform.resources,skewed.resources)

readr::write_delim(all.data,"./results/abundances_niche_model_complete.csv",delim = ";")
