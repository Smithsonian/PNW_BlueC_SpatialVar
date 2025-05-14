# Jags code to model spatial scaling in blue carbon soils

# Import packages
library(rjags)
library(tidyverse)
library(ggmcmc)
library(stringr)

# Import dataset
blue_c_data <- read_csv("data/NE-Pacific-NAT-stocks.csv")

# Prep dataset
names(blue_c_data)

blue_c_data_long <- blue_c_data %>% 
  mutate(core_n = 1:n()) %>% 
  # Decide what level variables we are going to include
  # Lv1 Ecoregion
  # KG climate zone
  # Estuary
  select(core_n, Ecosystem, Stk30, Stk50, Stk100, Lvl1EcoReg, KGzone, Estuary) %>% 
  # Pivot
  gather(key = "depth", value = "stock", -c(core_n, Ecosystem, KGzone, Lvl1EcoReg,Estuary)) %>% 
  filter(!is.na(stock)) %>% 
  arrange(core_n) %>% 
  # Log transform
  mutate(ln_stock = log(stock)) %>% 
  mutate(depth = as.numeric(str_remove_all(depth, "Stk"))) %>% 
  arrange(Ecosystem, depth)

blue_c_data_summary <- blue_c_data_long %>% 
  group_by(Ecosystem, depth) %>% 
  summarise(mean = mean(stock))

ggplot(blue_c_data_long, aes(x = stock)) +
  geom_density() +
  facet_grid(Ecosystem~depth, scale = "free")

ggsave("figs/supplemental/All_Stock_Distributions.jpg", height=8.5, width = 11)

ggplot(blue_c_data_long, aes(x = stock)) +
  geom_density(aes(fill = Estuary, color = Estuary),  alpha = 0.5) +
  facet_wrap(Ecosystem~depth, scale = "free") +
  theme(legend.position = "none") +
  ggtitle("Estuary")

ggsave("figs/supplemental/Estuary_Distributions.jpg", height=8.5, width = 11)

ggplot(blue_c_data_long, aes(x = stock)) +
  geom_density(aes(fill = KGzone, color = KGzone),  alpha = 0.5) +
  facet_grid(Ecosystem~depth, scale = "free") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("KG Zone")

ggsave("figs/supplemental/KG_Distributions.jpg", height=8.5, width = 11)

ggplot(blue_c_data_long, aes(x = stock)) +
  geom_density(aes(fill = as.character(Lvl1EcoReg), color = as.character(Lvl1EcoReg)),  alpha = 0.5) +
  facet_grid(Ecosystem~depth, scale = "free") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("Level-1 Ecoregion")

ggsave("figs/supplemental/LV1_Distributions.jpg", height=8.5, width = 11)

# Unique numeric id's
core_depth_table <- blue_c_data_long %>% 
  select(Ecosystem, depth) %>% 
  distinct_all() %>% 
  arrange(Ecosystem, depth) %>% 
  mutate(ecosystem_depth = 1:n())

# Unique numeric id's for tier 1, kg climate zone, estuary, habitat + depth
blue_c_data_jaggified <- blue_c_data_long %>% 
  left_join(core_depth_table)

# Check to see if all depth-eco's have at least 2 from each level of heirarchy
Estuary_check <- blue_c_data_jaggified %>% 
  select(ecosystem_depth, Estuary) %>% 
  distinct_all() %>% 
  group_by(ecosystem_depth) %>% 
  summarise(n = n())

(Estuary_check)

kg_check <- blue_c_data_jaggified %>% 
  select(ecosystem_depth, KGzone) %>% 
  distinct_all() %>% 
  group_by(ecosystem_depth) %>% 
  summarise(n = n())

(kg_check)

t1_check <- blue_c_data_jaggified %>% 
  select(ecosystem_depth, Lvl1EcoReg) %>% 
  distinct_all() %>% 
  group_by(ecosystem_depth) %>% 
  summarise(n = n())

(t1_check)

# Convert to numbers
blue_c_data_jaggified$Estuary <- as.numeric(factor(blue_c_data_jaggified$Estuary))
blue_c_data_jaggified$KGzone <- as.numeric(factor(blue_c_data_jaggified$KGzone))
blue_c_data_jaggified$Lvl1EcoReg <- as.numeric(factor(blue_c_data_jaggified$Lvl1EcoReg))

# Jags code

spatial_model <- "
model{

  # Set up uninformed priors
  # Iterate through all 15 depth and ecosystems
  for(m in 1:15) {
  
    # Priors for means each Ecosystem and depth combo
    mu[m] ~ dunif(0, 500)
    
    # Priors for each level of tau
    # tier 1 climate zone - for each ecosystem x depth
    tau_tier1[m] ~ dgamma(0.001, 0.001)
    
    # Covert precisions to variances
    var_tier1[m] <- 1/tau_tier1[m]
    
    # kg climate zone - for each ecosystem x depth
    tau_kg[m] ~ dgamma(0.001, 0.001)
    var_kg[m] <- 1/tau_kg[m]
    
    # estuary for each ecosystem x depth
    tau_est[m] ~ dgamma(0.001, 0.001)
    var_est[m] <- 1/tau_est[m]
    
    # core-level for each ecosystem x depth
    tau_core[m] ~ dgamma(0.001, 0.001)
    var_core[m] <- 1/tau_core[m]
    
    # tier 1 climate zones l
    # for l in unique T1 climate zone, ecosystem, depth combo
    # start l
    for (l in 1:L) {
        beta_tier1[l,m] ~ dnorm(0, tau_tier1[m]) 
    } # end l
    
    # kg climate zones k
    # for iterator in unique KG climate zone, ecosystem, depth combo
    # start k
    for(k in 1:K) {
      beta_kg[k, m] ~ dnorm(0, tau_kg[m]) 
    } # end k
    
    # estuary j
    # for iterator in unique estuary, ecosystem-depth combo
    # start j
    for (j in 1:J) {
      beta_est[j,m] ~ dnorm(0, tau_est[m]) 
    }
    # end j
    
  
  }  # end m
  
  # Process model
  for (i in 1:I) {
    y[i] ~ dnorm(mu[ecosystem_depth[i]] + beta_tier1[tier1[i], ecosystem_depth[i]] +  beta_kg[kg[i], ecosystem_depth[i]] +  beta_est[est[i], ecosystem_depth[i]], tau_core[ecosystem_depth[i]]) T(0,)
  }

}
"
# Run jags code
data_list <- list(y = blue_c_data_jaggified$stock,
                  ecosystem_depth = blue_c_data_jaggified$ecosystem_depth,
                  tier1 = blue_c_data_jaggified$Lvl1EcoReg,
                  kg = blue_c_data_jaggified$KGzone,
                  est = blue_c_data_jaggified$Estuary,
                  I = nrow(blue_c_data_jaggified),
                  J = max(blue_c_data_jaggified$Estuary),
                  K = max(blue_c_data_jaggified$KGzone),
                  L = max(blue_c_data_jaggified$Lvl1EcoReg))

j.model   <- jags.model(file = textConnection(spatial_model),
                        data = data_list,
                        n.chains = 4,
                        n.adapt=2000,
                        inits = list(list(mu = blue_c_data_summary$mean + rnorm(15, 0, 1)),
                                     list(mu = blue_c_data_summary$mean + rnorm(15, 0, 1)),
                                     list(mu = blue_c_data_summary$mean + rnorm(15, 0, 1)),
                                     list(mu = blue_c_data_summary$mean + rnorm(15, 0, 1))))

var.out   <- coda.samples(model = j.model,
                          variable.names = c("mu",
                                             "var_core",
                                             "var_est",
                                             "var_kg",
                                             "var_tier1"),
                          , n.iter = 4000)


# Make tidy jags
tidyJags<- ggs(var.out) 

# Check convergence
ggmcmc(tidyJags,
       plot = c("ggs_traceplot", "ggs_density"))

# Calculate proportion of cumulative variance
proportional_variance <- tidyJags %>% 
  filter(grepl("var", Parameter)) %>% 
  separate(Parameter, into = c("parameter", "ecosystem_depth"), sep = "\\[") %>% 
  mutate(ecosystem_depth = str_remove_all(ecosystem_depth, "\\]")) %>% 
  group_by(Chain, Iteration, ecosystem_depth) %>% 
  mutate(cumulative_var = sum(value)) %>% 
  ungroup() %>% 
  mutate(proportional_variance = value / cumulative_var * 100)

# Summarise as mean, median, se
summarized_var <- proportional_variance %>% 
  group_by(parameter, ecosystem_depth) %>% 
  summarise(mean = mean(proportional_variance),
            median = median(proportional_variance),
            sd = sd(proportional_variance),
            upper_CI = quantile(proportional_variance, 0.84),
            lower_CI = quantile(proportional_variance, 0.16))

# Graph
graph_var <- summarized_var %>% 
  mutate(ecosystem_depth = as.numeric(ecosystem_depth)) %>% 
  left_join(core_depth_table) %>% 
  mutate(Ecosystem = recode(Ecosystem, "EM"="Emergent marsh",
                            "MG" = "Mangrove",
                            "FL" = "Tidal flat",
                            "SG" = "Seagrass",
                            "TS" = "Tidal swamp"),
         Ecosystem = factor(Ecosystem, levels = c("Tidal flat", 
                                                  "Seagrass",
                                                  "Emergent marsh",
                                                  "Mangrove",
                                                  "Tidal swamp")),
         parameter = str_remove_all(parameter, "var_"),
         depth = as.character(depth),
         depth = recode(depth, 
                        "30" = "0-30 cm",
                        "50" = "0-50 cm",
                        "100" = "0-100 cm"),
         depth = factor(depth, levels = c("0-30 cm",
                                          "0-50 cm",
                                          "0-100 cm")),
         parameter = recode(parameter, 
                            "core"="within estuary",
                            "est" = "estuary",
                            "kg"= "Köppen–Geiger climate zone",
                            "tier1"="Level-1 ecoregion"),
         parameter = factor(parameter, levels = c("within estuary",
                                                   "estuary",
                                                   "Köppen–Geiger climate zone",
                                                   "Level-1 ecoregion")))
  
write_csv(graph_var, "summary_table.csv")
 
ggplot(graph_var, aes(x = parameter, y = median)) +
  geom_point() +
  geom_segment(aes(y = lower_CI, yend = upper_CI, xend = parameter)) +
  facet_grid(depth~Ecosystem) +
  ylab("Proportional Variance (%)") +
  xlab("Spatial Level (finer to coarser scale)") +
  # Chris's specs 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave("figs/main/Satial_variance.jpg", height=188*0.66, width = 188, unit = "mm")



