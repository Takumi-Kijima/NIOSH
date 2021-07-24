###########################################
#### NIOSH ten Berge Exponent Analysis ####
####   Takumi Kijima & Matthew Snyder  ####
####        STA 660, FA 2020           ####
###########################################

library(readr)
library(tidyverse)
library(ggfortify)
library(rstan)
library(bayesplot)
library(HDInterval)
library(tcltk)
library(gridExtra)
######################################################################################### Read in data
niosh0 <- read_csv("NIOSH.csv")

# Data manipulation and transformation
niosh <- niosh0 %>% 
  mutate(LC50n = as.numeric(LC50),
         Time_min = case_when(Time_Units == "hr" ~ as.numeric(Time)*60,
                              Time_Units == "min" ~ as.numeric(Time))) %>% 
  filter(!is.na(LC50n) & !is.na(Time_min)) %>% 
  mutate(logLC50 = log(LC50n),
         logTime = log(Time_min))

# Summary stats of working chemicals
summary <- niosh %>% 
  group_by(Chemical, Species) %>% 
  summarize(n=n(),
            unique_times=length(unique(Time_min)))

summary_table <- summary %>% 
  pivot_wider(values_from=c(n, unique_times), names_from=Species, id_cols=Chemical)

##################################################### Filter only those chemicals with 3 or more rat obs and 2 or more times
good_chem <- summary %>% 
  filter(Species=="Rat",
         n >= 3,
         unique_times >= 2)
good_niosh_2 <- niosh %>% 
  filter(Species == "Rat" & Chemical %in% good_chem$Chemical)

# Log transformed plot
ggplot(good_niosh_2) +
  geom_point(aes(x=logTime, y=logLC50, color=Chemical)) +
  geom_line(aes(x=logTime, y=logLC50, color=Chemical), alpha=.3) +
  facet_wrap(vars(Chemical)) +
  labs(x="Log of Time (minutes)", y="Log of LC50 (ppm)", title="LC50 Values for Rats over Time by Chemical") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5))

# Original data plot
ggplot(good_niosh_2) +
  geom_point(aes(x=Time_min, y=LC50)) +
  facet_wrap(vars(Chemical))


#################################################################################################### Multispecies
multiple <- summary %>% 
  filter(n >= 3, unique_times >=2) %>% 
  group_by(Chemical) %>% 
  filter(n() > 1)
mult_niosh <- niosh %>% 
  left_join(multiple, by=c("Species", "Chemical")) %>% 
  filter(!is.na(unique_times))


ggplot(mult_niosh) +
  geom_point(aes(x=logTime, y=logLC50, color=Species)) +
  geom_line(aes(x=logTime, y=logLC50, color=Species), alpha=.3) +
  facet_wrap(vars(Chemical), scales="free_x") +
  labs(x="Log of Time (minutes)", y="Log of LC50 (ppm)", title="LC50 Values over Time for Chemicals with Multiple Species") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5))


#################################################################################################### Filter chemicals with 3 rat obs from 3 different times
good_chem_3 <- summary %>% 
  filter(Species=="Rat",
         n >= 3,
         unique_times >= 3)
good_niosh <- niosh %>% 
  filter(Species == "Rat" & Chemical %in% good_chem_3$Chemical)

# Log transformed plot
ggplot(good_niosh) +
  geom_point(aes(x=logTime, y=logLC50, color=Chemical)) +
  geom_line(aes(x=logTime, y=logLC50, color=Chemical), alpha=.3) +
  facet_wrap(vars(Chemical)) +
  labs(x="Log of Time (minutes)", y="Log of LC50 (ppm)", title="LC50 Values for Rats over Time by Chemical") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5))



#################################################################################################### Regression
# Function to fit loglinear model
loglm <- function(chem) {
  data <- good_niosh %>% 
    filter(Chemical == chem)
  
  mod <- lm(logLC50 ~ logTime, data=data)
  mod
}

# Fit models to all chemicals in  good_niosh
models <- list()
for (chem in good_niosh$Chemical) {
  chem_name <- paste0(chem)
  models[[chem_name]] <- loglm(chem)
}

# Residuals of models
residuals <- list()
for (i in 1:length(models)) {
  chem_name <- names(models)[i]
  residuals[[chem_name]] <- autoplot(models[[i]])
}

residuals[[1]]
residuals[[2]]
residuals[[3]]
residuals[[4]]
residuals[[5]]  # violations
residuals[[6]]
residuals[[7]] # CV violation
residuals[[8]]
# residuals[[9]]
# residuals[[10]]
# residuals[[11]]
# residuals[[12]]
# residuals[[13]]


# beta estimates and n calculations
coeff <- data.frame("Chemical"=character(), "beta1"=numeric(), "sd"=numeric(), 
                    "l_beta"=numeric(), "u_beta"=numeric(), 
                    "n"=numeric(), "l_n"=numeric(), "u_n"=numeric())
for (i in 1:length(models)) {
  coeff[i,1] <- names(models)[i]
  
  mod_sum <- summary(models[[i]])
  # Beta coefficient
  coeff[i,2] <-mod_sum$coefficients[2,1]
  # Beta sd
  coeff[i,3] <- mod_sum$coefficients[2,2]
  # Beta lower bound
  coeff[i,4] <- confint(models[[i]])[2,1]
  # Beta upper bound
  coeff[i,5] <- confint(models[[i]])[2,2]
  # n (ten Berge Exponent)
  coeff[i,6] <- -1/coeff[i,2]
}




#################################################################################################### Bootstrapping
### https://www.statmethods.net/advstats/bootstrapping.html
library(boot)
set.seed(123)

# Bootstapping function for CI of tenBerge exponent
tenBerge <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  mod <- lm(formula, data=d)
  
  beta1 <-mod$coefficients[2]
  
  return(-1/beta1)
}


for (i in 1:length(models)) {  
  data <- good_niosh %>% 
    filter(Chemical == names(models)[i])
  
  results <- boot(data=data, statistic=tenBerge,
                  R=1000, formula=logLC50 ~ logTime)
  
  # view results
  #results
  #plot(results)
  
  # get 95% confidence interval
  ci <- boot.ci(results, type="bca")
  coeff[i,7] <- ci[[4]][4]
  coeff[i,8] <- ci[[4]][5]
}


ggplot(coeff, aes(x=Chemical, color=Chemical)) +
  geom_errorbar(aes(ymin = l_n, ymax = u_n), width = 0.2) +
  geom_point(aes(y=n)) +
  geom_hline(aes(yintercept=.85)) +
  geom_hline(aes(yintercept=3)) +
  theme_bw()








