# Dose-response models code #########################################

# The aim of this script is to compare the predicted responses for dose-response (at discharge time point) vs. time-course models at TREATEMENT-LEVEL

# Libraries required
library(readxl)
library(MBNMAdose)
library(tidyverse)
library(ggthemes)
library(patchwork)
library(cowplot)
library(igraph)
library(esc)
library(purrr)
library(viridis) 


# Load dataset
data <- read_xlsx(path = paste0(getwd(), "/hosp_dataset.xlsx"), # name of the dataset file 
                  col_names = TRUE, 
                  sheet = "hospital_data_agent_MBNMAdose") # name of the sheet

glimpse(data)

# Effect sizes (Standardised Mean Differences) estimation from means and standard deviations: change from baseline approach
data_es <- do.call(rbind, lapply(1:nrow(data), function(i) {
  esc_mean_sd(
    grp1m = data$Me[i], 
    grp1sd = data$SDe[i], 
    grp1n = data$Ne[i],
    grp2m = data$Mc[i], 
    grp2sd = data$SDc[i], 
    grp2n = data$Nc[i],
    study = data$studyID[i], 
    es.type = "g"
  )
})) %>% as.data.frame()

# Data cleaning
data <- data %>% 
  select(studyID, Ne, dose, agent, lower_better, Me, SDe) %>% 
  rename(n = Ne) %>% 
  mutate(y = unlist(data_es$es), 
         se = unlist(data_es$se))


# Network creation
net <- mbnma.network(data)

plot(net, v.color = "agent", label.distance = 3, remove.loops = TRUE) # treatment-level

plot(net, level = "agent", label.distance = 3.5) # agent-level

summary(net) # network description


# Examining dose-response relationships
## First and foremost, it is very useful to visualize the effects estimated using a standard (i.e., “split”) Network Meta-Analysis to explore the potential dose-response functions that could work to the responses pattern
nma.data <- nma.run(net, method = "random")

nma.data # results of the NMA

plot(nma.data) # plot of NMA estimates plotted by dose 

### There appears to be a dose-response relationship and it also appears to be non-linear


## We can also conduct an Unrelated Mean Effects (UME) model to compare with consistent model fit and see whether global consistency assumption holds
nma.ume <- nma.run(net, method = "random", UME = TRUE)

data.frame(deviance = c(nma.data$jagsresult$BUGSoutput$summary[12, 1], 
                        nma.ume$jagsresult$BUGSoutput$summary[67, 1]), 
           SD = c(nma.data$jagsresult$BUGSoutput$summary[13, 1], 
                  nma.ume$jagsresult$BUGSoutput$summary[68, 1])) # if values are similar, consistency holds


# Analysis using mbnma.run to perform dose-response candidate functions

## Emax dose-response function using common-treatment effects
emax <- mbnma.run(net, fun=demax(emax="rel", ed50="rel"),
                  method="common")

## Quadratic dose-response function using common-treatment effects
func <- ~ (beta.1 * dose) + (beta.2 * (dose^2))
quad <- mbnma.run(net, fun=duser(fun=func), method="common")

## Restricted cubic spline dose-response function using common-treatment effects
knots <- c(0.1, 0.5, 0.9)
rcs <- mbnma.run(net, fun = dspline(type = "bs", knots = knots), 
                 method = "common")

### Fit comparison
data.frame(deviance = c(emax$BUGSoutput$summary[7, 1],
                        quad$BUGSoutput$summary[7, 1],
                        rcs$BUGSoutput$summary[7, 1]),
           DIC = c(emax$BUGSoutput$DIC,
                   quad$BUGSoutput$DIC,
                   rcs$BUGSoutput$DIC),
           row.names = c("Emax", "Quadratic", "rcs")) # quadratic and rcs dose-response functions presented similar fits

## Quadratic dose-response function using random-treatment effects
random_q <- mbnma.run(net, fun=duser(fun=func), method="random")

## Restricted cubic spline dose-response function using random-treatment effects
random_rcs <- mbnma.run(net, fun = dspline(type = "bs", knots = knots), 
                        method = "random")

## Fit comparison
data.frame(deviance = c(quad$BUGSoutput$summary[7, 1],
                        rcs$BUGSoutput$summary[7, 1],
                        random_q$BUGSoutput$summary[7, 1],
                        random_rcs$BUGSoutput$summary[7, 1]),
           DIC = c(quad$BUGSoutput$DIC,
                   rcs$BUGSoutput$DIC,
                   random_q$BUGSoutput$DIC,
                   random_rcs$BUGSoutput$DIC),
           SD = c(NA,
                  NA,
                  random_q$BUGSoutput$summary[8, 1],
                  random_rcs$BUGSoutput$summary[8, 1]),
           row.names = c("Quadratic common", "rcs common", 
                         "Quadratic random", "rcs random")) # common-treatment effects had better fit 


# Post-estimations

## Deviance plots 
devplot(quad, plot.type = "box") # a point at dose 100 in Multicomponent intervention seems to not fit 

## Fitted values plot
fitplot(quad)


## Forest plot
plot(quad) # beta coefficients’ distribution
plot(rcs) # as sensitivity analysis


## Predicted responses 
### Specify the doses
doses <- list("Ambulation" = seq(0, 250), 
              "Mobility" = seq(0, 250),
              "Multicomponent" = seq(0, 250))

### Predict exact doses and estimates assuming E0 response = 0
pred <- predict(quad, E0 = 0, exact.doses = doses)

plot(pred) # predicted responses visualisation



#### Sensitivity analysis predicting responses using the rcs model
### Considering that Emax also had a good model fit, we can visualize the predicted responses using this function
plot(predict(emax, E0 = 0, exact.doses = doses))
## These responses make biological sense, where patients will benefit from physical activity up to a point in which the extra benefits of moving more are not significant

# see responses
data.frame(summary(emax)) %>%
  filter(dose %in% c(0, 50, 100, 150, 200, 250))

## Ranking of effectiveness
### First, we have to conduct a model with the predicted responses for specified doses at 50, 100, 150, 200 and 250 METs-min/day
pred_rank <- predict(quad, E0 = 0, 
                     exact.doses = list("Ambulation" = c(0, 50, 100, 150, 200, 250),
                                        "Mobility" = c(0, 50, 100, 150, 200, 250),
                                        "Multicomponent" = c(0, 50, 100, 150, 200, 250)))
                     
### Conduct ranking
rank <- rank(pred_rank, direction = 1)
                     
plot(rank) # plotted ranking probabilities
                
   
                    
# Visualisation of results
## Tile plots

# Ambulation plot
m1 <- data.frame(summary(pred)) %>% 
  filter(agent %in% c("Ambulation")) %>%
  ggplot(aes(x = dose, y = agent, fill = mean)) +
  # the color indicates the effect, corresponding from purple (no effect) to yellow (the greatest effect)
  viridis::scale_fill_viridis(discrete = FALSE) +
  geom_tile() +
  # PA dose where a minimal effect is detected for the patients
  geom_vline(xintercept = 74, col = "red", size = 1.1, alpha = 0.7) +
  labs(y = "", x = "", fill = expression(gamma), 
                            title = "Ambulation") + 
  theme_minimal_vgrid() +
  theme(axis.text.y = element_blank(),
                             legend.text = element_text(size = 9),
                             legend.box.spacing = unit(0.3, "cm"),
                             legend.key.size = unit(0.4, "cm")) +
  scale_x_continuous(breaks = seq(0, 250, 25))


# Mobility plot
m2 <- data.frame(summary(pred)) %>% 
  filter(agent %in% c("Mobility")) %>%
  ggplot(aes(x = dose, y = agent, fill = mean)) +
  viridis::scale_fill_viridis(discrete = FALSE) +
  geom_tile() +
  labs(y = "", x = "", fill = expression(gamma),
                            title = "Mobility") + 
  theme_minimal_vgrid() +
  theme(axis.text.y = element_blank(),
                             legend.text = element_text(size = 9),
                             legend.box.spacing = unit(0.5, "cm"),
                             legend.key.size = unit(0.4, "cm")) +
  scale_x_continuous(breaks = seq(0, 250, 25))
           

# Multicomponent plot          
m3 <- data.frame(summary(pred)) %>% 
  filter(agent %in% c("Multicomponent")) %>%
  ggplot(aes(x = dose, y = agent, fill = mean)) +
  viridis::scale_fill_viridis(discrete = FALSE) +
  geom_tile() +
  geom_vline(xintercept = 226, col = "red", alpha = 0.8, size = 1.1) +
  labs(y = "", x = "Exercise dose (METs-min/day)", 
       fill = expression(gamma),
                            title = "Multicomponent") + 
  theme_minimal_vgrid() +
  theme(axis.text.y = element_blank(),
                             legend.text = element_text(size = 9),
                             legend.box.spacing = unit(0.5, "cm"),
                             legend.key.size = unit(0.4, "cm")) +
  scale_x_continuous(breaks = seq(0, 250, 25))

# plot!                     
plot_grid(m1, m2, m3, ncol = 1) # combine plots

