# IEP Workshop Training: Simple linear regression
# Author: Braden Elliott
# 23 March 2022

# RESEARCH QUESTION: Are longer fish heavier? Specifically, can we predict a fish's weight from its length?

#####################
# Prepare R session #
#####################

# install.packages(c("FSA","FSAdata", "ggplot2")) # Only uncomment and run this line if packages not already installed.
#library(FSA) # library([name]) loads the named package into the current R session, so functions it supports may be used.
library(FSAdata)
library(ggplot2)

###############################
# Import and and prepare data #
###############################

data(RuffeSLRH92) # read in data on ruffe, a species of fish
str(RuffeSLRH92) # look at the structure of the data
summary(RuffeSLRH92) # look at quantitative summaries of the data
head(RuffeSLRH92)
ruffe2 <- subset(RuffeSLRH92,!is.na(weight) & !is.na(length)) # drop lines of the dataset which lack one or both of the variables of interest
plot(RuffeSLRH92$length, RuffeSLRH92$weight) # not linear, need to transform to make linear.
ruffe2$logL <- log(ruffe2$length) # log transform ruffe length
ruffe2$logW <- log(ruffe2$weight) # log transform ruffe weight
plot(ruffe2$logL,ruffe2$logW) # visualize relationship between log(length) and log(weight) of ruffe
ruffe3 <- subset(ruffe2,logW >= -0.5) # drop lines of the dataset under 0.5 log(weight) since those data points are relatively unreliable for this analysis

################
# Analyze data #
################

ruffeLW <- lm(logW~logL,data=ruffe3) # run linear model, predicting log(weight) from log(length).

summary(ruffeLW) # show the summarized output of the linear model: p-values

##################
# Validate model #
##################

# Are the residuals (variation between prediction and actual value for each data point) normally distributed?
hist(residuals(ruffeLW)) # rough visualization via histogram (looking for bell shape)
qqnorm(residuals(ruffeLW)) # precise visualization via quantile-quantile plot (looking for straight line from origin)
plot(ruffe3$logL,residuals(ruffeLW)) # check for funneling or any other shapes, which would indicate heteroscedasticity
plot(ruffe3$length,residuals(lm(weight~length,data=ruffe3))) # example of heteroscedasticity in untransformed ruffe data -- this is why we log-transformed during data prep above.

#################
# Present model #
#################

ggplot(ruffe3,aes(logL, logW)) + # initiate a plot command using ggplot2, based on the ruffe3 dataset and defining log(length) as the X (horizontal) axis since it is the predictor and log(weight) as the y (vertical) axis since it is the response.
  geom_point() + # add points to the plot, which will read by default from the variables used to define the axes.
  geom_smooth(method='lm', se=FALSE, color='turquoise4') + # add linear model trendline to the plot, suppressing the standard error margin
  theme_minimal() + # call in an aesthetically-pleasing built-in design theme for rending the plot
  labs(x='log Length', y='log Weight', title='Relationship between Length and Weight') + # Add axis labels and a title
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) # define a couple details for rendering text, and by not ending the line with another "+" the ggplot command is closed.

################
# Bonus tracks #
################

# what if multiple predictors affect the response?
ruffeALW <- lm(logW~logL+age,data=ruffe3) # run linear model, predicting log(weight) from log(length) AND age
summary(ruffeALW) # in this case, age does not add any useful information to the model (small effect, low significance).Note data points dropped due to missingness, meaning some lines were missing age data.

# Can I make predictions from the model on new data?
# Yes. The "predict" function is used for this purpose, use ?predict to learn about it.
# Remember to back-transform if using transformed variables in the model!

# Can I do anything else to validate my model?
# install.packages("DAAG") # install package if necessary
library(DAAG) # pull up a package to run k-fold validation, which subdivides your data randomly to see if each subset has the same pattern.
cv.lm(data=ruffe3, ruffeLW, m=10) # look at ten subsets ("folds") modeled independently and plotted together
