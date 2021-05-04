# Load required packages and data
library(ape)
library(brms)
phylo <- read.nexus("phylogeny.nex") 
data <- read.csv("data.csv", header=TRUE)

# Specify brms model
pf.mod <- bf(pf ~ humidity + mass + richness + (1|gr(species, cov = A)))
mass.mod <- bf(mass ~ humidity + (1|gr(species, cov = A)))
richness.mod <- bf(richness ~ humidity + (1|gr(species, cov = A)))
priors <- c(set_prior("normal(0,2)", class = "b", coef = "mass", resp = "pf"), 
            set_prior("normal(0,2)", class = "b", coef = "richness", resp = "pf"),
            set_prior("normal(0,2)", class = "b", coef = "humidity", resp = "pf"),
            set_prior("normal(0,2)", class = "b", coef = "humidity", resp = "mass"),
            set_prior("normal(0,2)", class = "b", coef = "humidity", resp = "richness"))

# Prepare phylogenetic correlation matrix
rownames(data) <- data$species
humidity <- drop.tip(humidity, setdiff(humidity$tip.label, data$species))
data <- data[order(match(data$species,humidity$tip.label)),]
A <- ape::vcv.phylo(humidity, corr = TRUE)

# Run brms model
fit <- brm(pf.mod + mass.mod + richness.mod + set_rescor(FALSE), 
    data=data, data2 = list(A = A), sample_prior = TRUE,
    cores=2, chains = 2, iter = 10000, warmup = 2000, thin = 4,
    control = list(adapt_delta = 0.99), prior = priors)

# Extract posterior distributions for parameters of interest
parameters <- c("b_pf_humidity", "b_mass_humidity", "b_pf_mass", "b_richness_humidity","b_pf_richness", "sd_species__pf_Intercept", "sd_species__mass_Intercept", "sd_species__richness_Intercept", "sigma_pf", "sigma_mass", "sigma_richness")
posterior <- posterior_samples(fit, pars = parameters)

# Calculate posterior distributions for indirect effects
indirect_mass <- posterior$b_mass_humidity * posterior$b_pf_mass
indirect_richness <- posterior$b_richness_humidity * posterior$b_pf_richness


# Calculate posterior distributions for phylogenetic signal
lambda_pf <- parameters$sd_species__pf_Intercept^2 * (parameters$sd_species__pf_Intercept^2 + parameters$sigma_pf)
lambda_mass <- parameters$sd_species__mass_Intercept^2 * (parameters$sd_species__mass_Intercept^2 + parameters$sigma_mass)
lambda_richness <- parameters$sd_species__richness_Intercept^2 * (parameters$sd_species__richness_Intercept^2 + parameters$sigma_richness)