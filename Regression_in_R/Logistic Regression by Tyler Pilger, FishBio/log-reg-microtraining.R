##------------------------------------------------------------------------ BEGIN
# File-Name:        log-reg-microtraining.r
# Version:          v.1
# Date:             2022-3-1
# Authors:          Tyler Pilger, Matt Peterson
# Email(s):         tylerpilger@fishbio.com, mattpeterson@fishbio.com 
# Description:      Companion script for the 2022 IEP Workshop statistics micro-training session on
#                   logistic regression, how to model two category, binary, or proportional data.
#                   

# Load libraries
library(dplyr)
library(ggplot2)

### Import Data
# Splittail (Pogonichthys macrolepidotus) records from the Delta Juvenile Fish Monitoring Program (USFWS)
# at all stations during the years 2001 to 2010.
splt <- readRDS(file = 'DJFMP_seine_2001-2010_splt.rds')

# Filter to include months of Feb-June
splt <- filter(splt, Month %in% c(2:6))

# Plot the numbers of splittail caught against water temperature
ggplot(splt, aes(Wtemp, Count)) +
  geom_point()

# Create vector for presence/absence
splt$PA <- ifelse(splt$Count>0, 1, 0)

# Plot presence/absence data
ggplot(splt, aes(Wtemp, PA)) +
  geom_point() +
  labs(x='Water Temperature (C)', y='Presence-Absence')

# Create new dataset with proportion of stations present during each month of each year and the mean
# water temperature measured across stations
splt.prop <- group_by(splt, Year, Month) %>%
  summarise(., n = length(unique(StationCode)),
            nPres = length(unique(StationCode[which(PA==1)])),
            mnTemp = mean(Wtemp, na.rm=T)) %>%
  mutate(., propPres=nPres/n)

# Plot the proportion of sites that splittail were present against mean water temperature
ggplot(splt.prop, aes(mnTemp, propPres)) +
  geom_point() +
  labs(x='Mean Water Temperature (C)', y='Proportion of Stations with Splittail')


###-- Modeling presence/absence as a function of temperature

mod1 <- glm( PA ~ Wtemp, data = splt, family = 'binomial')
summary(mod1)
range(splt$Wtemp, na.rm=T)
newdat <- data.frame(Wtemp = seq(7, 28, length.out=50))

PredPA <- predict(mod1, newdata = newdat, type = 'response', se.fit = T)

newdat$PredPA <- PredPA$fit
newdat$PredPA.se <- PredPA$se.fit

# Plot predicted probability of splittail being present at a site as a function of water temperature 
ggplot(splt, aes(Wtemp, PA)) + 
  geom_point() +
  geom_line(data = newdat, aes(Wtemp, PredPA), size = 1.5) +
  geom_line(data = newdat, aes(Wtemp, PredPA + PredPA.se*1.96), size = .75, color='grey') +
  geom_line(data = newdat, aes(Wtemp, PredPA - PredPA.se*1.96), size = .75, color='grey') +
  labs(x='Water Temperature (C)', y='P(Presence)')

# Using plot() with presence-absence data is tough to interpret
op <- par(mfrow=c(2,2))
plot(mod1, ask = F)
par(op)

# Graphical test for fit of the model
cutT <- cut(splt$Wtemp, 5)
tapply(splt$PA, cutT, sum)
table(cutT)
# empirical probabilities of the means
probs <- tapply(splt$PA, cutT, sum) / table(cutT)
tempMeans <- tapply(splt$Wtemp, cutT, mean)

# Calculate standard error of the binomial proportions
se <- sqrt(probs * (1-probs)/table(cutT))

probMeans <- data.frame(Wtemp = as.vector(tempMeans),
                        probs = as.vector(probs),
                        se = as.vector(se))

ggplot(splt, aes(Wtemp, PA)) + 
  geom_rug(data = splt[splt$PA==0,], sides = 'b') +
  geom_rug(data = splt[splt$PA==1,], sides = 't') +
  geom_line(data = newdat, aes(Wtemp, PredPA), size = 1.5) +
  geom_line(data = newdat, aes(Wtemp, PredPA + PredPA.se*1.96), size = .75, color='grey') +
  geom_line(data = newdat, aes(Wtemp, PredPA - PredPA.se*1.96), size = .75, color='grey') +
  geom_errorbar(data = probMeans, aes(Wtemp, ymin=probs-se, ymax=probs+se), width=.2, inherit.aes = F) +
  geom_point(data = probMeans, aes(Wtemp, probs), size = 2, inherit.aes = F) +
  labs(x='Water Temperature (C)', y='P(Presence)')

###-- Modeling proportion of sites as a function of temperature
# response variable is proportion of sites where present, 
mod2 <- glm( propPres ~ mnTemp, data = splt.prop, family = 'binomial', weights = n)
summary(mod2)

# response variable is successes (present) and failures absent
splt.prop$nAbs <- splt.prop$n-splt.prop$nPres

mod3 <- glm(cbind(nPres,nAbs) ~ mnTemp, data = splt.prop, family = 'binomial')
summary(mod3)
## NOTE: these show signs of overdispsersion, residual deviance >> residual degrees of freedom, could try quasibinom
mod4 <- glm(cbind(nPres,nAbs) ~ mnTemp, data = splt.prop, family = 'quasibinomial')
summary(mod4)

# Get predicted probabilities from the model
range(splt.prop$mnTemp) # find range of mean water temperature
newdat <- data.frame(mnTemp = seq(10, 22, length.out=50))

PredPA <- predict(mod4, newdata = newdat, type = 'response', se.fit = T)

newdat$PredPA <- PredPA$fit
newdat$PredPA.se <- PredPA$se.fit

# Plot predicted and observed
ggplot(splt.prop, aes(mnTemp, propPres)) +
  geom_point() +
  geom_line(data = newdat, aes(mnTemp, PredPA), size = 1.5) +
  geom_line(data = newdat, aes(mnTemp, PredPA + PredPA.se*1.96), size = .75, color='grey') +
  geom_line(data = newdat, aes(mnTemp, PredPA - PredPA.se*1.96), size = .75, color='grey') +
  labs(x='Mean Water Temperature (C)', y='Proportion of Stations with Splittail')


# Check model fit using plot function
op <- par(mfrow=c(2,2))
plot(mod4, ask = F)
par(op)


