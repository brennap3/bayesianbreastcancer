# install.packages("caret")
# install.packages("FactoMineR")
# install.packages("R2jags")

library(tidyverse)
require(FactoMineR)
require(ggplot2)
library(caret)
library(R2jags)
library(plotly)



#getwd()
setwd("C:/Users/admin/Documents/bayesianbreastcancer")
breast.cancer<-read.csv("breast-cancer.data",header=FALSE,sep=",")
head(breast.cancer)
summary(breast.cancer)
#
# 7. Attribute Information:
#   1. Class: no-recurrence-events, recurrence-events
# 2. age: 10-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80-89, 90-99.
# 3. menopause: lt40, ge40, premeno.
# 4. tumor-size: 0-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44,
# 45-49, 50-54, 55-59.
# 5. inv-nodes: 0-2, 3-5, 6-8, 9-11, 12-14, 15-17, 18-20, 21-23, 24-26,
# 27-29, 30-32, 33-35, 36-39.
# 6. node-caps: yes, no.
# 7. deg-malig: 1, 2, 3.
# 8. breast: left, right.
# 9. breast-quad: left-up, left-low, right-up,	right-low, central.
# 10. irradiat:	yes, no.

colnames(breast.cancer)<-c("Class","age","menopause","tumor.size","inv.nodes","node.caps",
                           "degree.malignant","breast","breast.quadrant","irradiation")


breast.cancer$degree.malignant<- breast.cancer$degree.malignant %>% as.character() %>% as.factor()
##
#### eexplore the data
##
# number of categories per variable
cats = apply(breast.cancer, 2, function(x) nlevels(as.factor(x)))
cats
mca1 = MCA(breast.cancer, graph = FALSE)
mca1_vars_df = data.frame(mca1$var$coord, Variable = rep(names(cats),
                                                         cats))
mca1_obs_df = data.frame(mca1$ind$coord)
# plot of variable categories
ggplot(data = mca1_vars_df, aes(x = Dim.1, y = Dim.2, label = rownames(mca1_vars_df))) +
  geom_hline(yintercept = 0, colour = "gray70") + geom_vline(xintercept = 0,
                                                             colour = "gray70") + geom_text(aes(colour = Variable)) + ggtitle("MCA plot of variables using R package FactoMineR")
##
#####
##
# contingency table
table(breast.cancer$irradiation,breast.cancer$Class)
# Stacked + percent
ggplot(breast.cancer %>% group_by(irradiation,Class) %>% tally(), aes(fill=irradiation, y=n, x=Class)) +
  geom_bar(position="fill", stat="identity") +ggtitle("stacked bar chart percent")
# as a stacked bar
ggplot(breast.cancer %>% group_by(irradiation,Class) %>% tally(), aes(fill=irradiation, y=n, x=Class)) +
  geom_bar(position="stack", stat="identity") +ggtitle("stacked bar chart")
# there appears to be a strong association between irradiation and recurrence events as oppossed
# to a stron association with  no irradiation and no-recurrence events
#node-cap-no and Class
table(breast.cancer$node.caps,breast.cancer$Class)
ggplot(breast.cancer %>% group_by(node.caps,Class) %>% tally(), aes(fill=node.caps, y=n, x=Class)) +
  geom_bar(position="fill", stat="identity") +ggtitle("stacked bar chart percent")
ggplot(breast.cancer %>% group_by(node.caps,Class) %>% tally(), aes(fill=node.caps, y=n, x=Class)) +
  geom_bar(position="stack", stat="identity") +ggtitle("stacked bar chart")
#there is a stronger association with node cpass and recurrence events
# table(breast.cancer$irradiation,breast.cancer$Class)
#
#     no-recurrence-events recurrence-events
# no                   164                54
# yes                   37                31
modelString="
model {
   #Likelihood
  for (i in 1:nGroups) {
    obs[i] ~ dbin(p[i],n[i])
    p[i] ~ dbeta(a[i],b[i])
    a[i] ~ dgamma(1,0.01)
    b[i] ~ dgamma(1,0.01)
  }
}
"
#     no-recurrence-events recurrence-events
# no                   164                54
## > 54/(54+164)
## [1] 0.2477064
#The likelihood model indicates that the observed counts are modeled by a
# binomial distribution with a probability of p (fraction 54/(54+164)) from n trials (items= 218).
#The prior on each p is defined as a beta distribution with shape parameters a and b
#The hyperpriors for each a and b are drawn from imprecise (vague, flat) gamma distributions.
#As input, JAGS will need to be supplied with:
# the observed data (obs)
# the total number of observed items (n)
# the number of classification groups (nGroups)
# This all needs to be contained within a list object.
#The observed item frequencies
obs <- c(54, 164)
data.list <- list(obs = obs, n = c(218, 218), nGroups = 2)
data.list
##the nodes (estimated parameters) to monitor (return samples for)
##the number of MCMC chains (3)
##the thinning factor (10)
##the number of MCMC iterations - determined by the number of samples to save,
##the rate of thinning and the number of chains
mod = jags.model(textConnection(modelString), data=data.list, inits=NULL, n.chains=3)
update(mod, 1e3)
mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combined chains
plot(mod_sim)
gelman.diag(mod_sim)
autocorr.diag(mod_sim)
effectiveSize(mod_sim)
summary(mod_sim)
# Conclusions: Rhat and n.eff are now much better for the probability parameters.
# The estimated fractions for  non irradiated for recurrence and non-recurrence are:
#
#
#   2.5%    25%    50%    75%  97.5%
# p[1] 0.1957 0.2300 0.2486 0.2685 0.3097
# p[2] 0.6908 0.7309 0.7505 0.7703 0.8050
#
#
# Collectively, the fractions of 1/2, 1/2  do not fall within these ranges.
## now lets do it for
# table(breast.cancer$irradiation,breast.cancer$Class)
#
#     no-recurrence-events recurrence-events
# no                   164                54
# yes                   37                31
tot_samples=164+54+37+31
modelString1="
model {
   #Likelihood
  for (i in 1:nGroups) {
    obs[i] ~ dbin(p[i],n[i])
    p[i] ~ dbeta(a[i],b[i])
    a[i] ~ dgamma(1,0.01)  ## inf scale 
    b[i] ~ dgamma(1,0.01)  ## 
  }
}
"
obs <- c(54, 164,37,31)
data.list <- list(obs = obs, n = c(tot_samples, tot_samples,tot_samples,tot_samples), nGroups = 4)
mod1 = jags.model(textConnection(modelString1), data=data.list, inits=NULL, n.chains=3)
update(mod1, 1e3)
mod_sim1 = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)
mod_csim1 = as.mcmc(do.call(rbind, mod_sim1)) # combined chains
plot(mod_sim1)
gelman.diag(mod_sim1)
autocorr.diag(mod_sim1)
effectiveSize(mod_sim1)
summary(mod_sim1)
# Collectively, the fractions of 1/4, 1/4,1/4,1/4  do not fall within these ranges.
# there appears to be a strong association between irradiation and recurrence of cancer