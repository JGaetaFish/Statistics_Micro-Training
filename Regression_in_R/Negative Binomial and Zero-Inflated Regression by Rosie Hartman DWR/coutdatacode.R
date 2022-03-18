#COunt data demo
#Rosemary Hartman
#IEP workshop 2022
#rosemary.hartman@water.ca.gov


library(tidyverse)
library(MASS)
library(DHARMa)
library(glmmTMB)
library(lmtest)
library(pscl)
library(lubridate)
library(visreg)
library(emmeans)
library(effects)

######################################################
#Part 1. 
#non zero-inflated data

FMWT =  read_csv("2019.08.08_countmodels presentation/FMWT 1967-2021 Catch Matrix_updated_tidy.csv")
WYs = read_csv("2019.08.08_countmodels presentation/yearassignments.csv")
#add the water year type
FMWT = left_join(FMWT, WYs)

#subset just the striped bass data and total up all the bass caught per station
SB = filter(FMWT, Species %in% c("Striped Bass age-0", "Striped Bass age-1",
                                 "Striped Bass age-1", "Striped Bass age-2")) %>%
  group_by(Year, SampleDate, SurveyNumber, StationCode, Volume, Yr_type) %>%
  summarize(Catch = sum(Catch, na.rm = T))

#now calculate total bass per survey
SBtot = group_by(SB, Year, SurveyNumber, Yr_type) %>%
  summarize(samples = n(), Volume = sum(Volume, na.rm = T), Catch = sum(Catch, na.rm = T))

#quick plot
ggplot(SBtot, aes(x = Year, y = Catch, color = Yr_type))+ geom_point()

#quick histogram

ggplot(SBtot, aes(x = Catch)) + geom_histogram()

#What does log-transformation look like?
ggplot(SBtot, aes(x = log(Catch))) + geom_histogram()


#try a linear model on log-transformed data
sblm = lm(log(Catch+1)~  Year + Yr_type, data = SBtot) 
summary(sblm)
plot(sblm)

#GLM using a poisson distribution
lm2 = glm(Catch ~ Year + Yr_type, data = SBtot, family = "poisson")
summary(lm2)
plot(lm2)

#different type of diganostic plot
simOut <- simulateResiduals(fittedModel = lm2, n_sim = 200)
plot(simOut, asFactor = F)

#test for overdispersion
testDispersion(simOut)

#now let's do a negative binomial distribution

sbnb = glm.nb(Catch~  Year + Yr_type,  data = SBtot) 
summary(sbnb)

#we can also use the DHARMa package
simOut <- simulateResiduals(fittedModel = sbnb, n_sim = 200)
plot(simOut, asFactor = F)
testDispersion(simOut)

#look at effects of year type on the log scale
emm = emmeans(sbnb, specs = "Yr_type")

emm

#now on the response scale
emm2 = emmeans(sbnb, specs = "Yr_type", type = "response")

emm2

#pairwise comparisons between year types

summary(contrast(emm2, method="pairwise", adjust="none", type="response"), infer=c(TRUE, TRUE))
visreg(sbnb)
visreg(sbnb, scale = "response")

#using offsets
SBtot$CPUE = SBtot$Catch/SBtot$Volume*10000
SBtot2 = dplyr::filter(SBtot, !is.na(Volume), Volume != 0)

sbnb2 = glm.nb(CPUE ~  Year,  data = SBtot2) 
#NOPE

sbnb3 = glm.nb(Catch ~ Year +offset(log(SBtot2$Volume)), data = SBtot2)
summary(sbnb3)



#######################################################################
#Part 2 - Zero Inflated data

#This time lets look at the catch of Dleta smelt in each individual trawl
#instead of the summary from all teh stations per survey


FMWT_DS = dplyr::filter(FMWT, Species == "Delta Smelt")

#just look at data from 1990-2020 in Montezuma Slough
FMWT_DS2 = dplyr::filter(FMWT_DS, Year <= 2020, Year > 1990, StationCode %in% c(605, 606, 608)) %>%
  mutate(Cond = scale(ConductivityTop), Sec = scale(Secchi), Year2 = scale(Year),
         Temp = scale(WaterTemperature))

#Some quick plots of the data

ggplot(FMWT_DS2, aes(x = Year, y = Catch, color = Yr_type)) +geom_point()

#now a histogram 
ggplot(FMWT_DS2) + geom_histogram(aes(x = Catch))

#now a log-transformed histogram
#(Note that we have to add 1 to each value to deal with the zeros)

ggplot(FMWT_DS2) + geom_histogram(aes(x = log(Catch+1)))

#That is definitely not going to work on a normal distribution!

#Poisson distributions usually don't work for fish data, so I"m 
#going to skip straight to negative binomial to see if that will work

#scale the predictor variables
FMWT_DS2 = mutate(FMWT_DS2, Cond = scale(ConductivityTop), Sec = scale(Secchi), Year2 = scale(Year),
                  Temp = scale(WaterTemperature))

#Negative binomial distribution.
#I've added Station as a random factor. You will learn why this afternoon. 
dsnb2 = glmmTMB(Catch~  Cond + Sec + Year2+ (1|StationCode),  family = "nbinom1", data = FMWT_DS2) 
summary(dsnb2)

#check out diagnostic plots
simres = simulateResiduals(dsnb2)
plot(simres)

testDispersion(simres)
testZeroInflation(simres)

#Some red text on there. What do we do?
#We know Delta Smelt populations have really declined in recent years, so maybe
#all those extra zeros are correlated to year!

dsnb4 = glmmTMB(Catch~  Sec +Cond+ Temp + (1|StationCode), 
                zi =~  Year2,
                family = "nbinom1", data = FMWT_DS2, 
                control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))) 
summary(dsnb4)
dsres = simulateResiduals(dsnb4, plot = T)


testDispersion(dsres)
testZeroInflation(dsres)
#Better!





