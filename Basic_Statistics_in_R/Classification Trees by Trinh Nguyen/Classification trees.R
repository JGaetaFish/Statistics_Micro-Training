# The script is a companion to the classification tree session for the 2022 IEP
# workshop. The dataset in question will be the Chipps Island dataset. 
# We will be predicting the degree of catch of longfin smelt (LFS) using various
# environmental variables.

# Libraries needed
# Reading in data
library(readr)
# Data manipulations
library(dplyr)
# Graphics
library(ggplot2)
# Caret provides a nice workflow interface and supports many different models
library(caret)
# Decision tree
library(rpart)
# Alternative visualization of the decision tree
library(rpart.plot)

# Dataset: Chipps Island; Focus species: LFS ------------------------------

data <- read_csv("Chipps Island Trawls CHN & POD Species 2012-2021.csv",
                 col_types = cols(
                   Location = col_character(),
                   Station = col_character(),
                   Date = col_date(format = "%m/%d/%Y"),
                   Time = col_time(format = ""),
                   Method = col_character(),
                   TowNumber = col_double(),
                   TowDuration = col_double(),
                   TowDirection = col_character(),
                   Volume = col_double(),
                   Secchi = col_double(),
                   DO = col_double(),
                   WaterTemp = col_double(),
                   Turbidity = col_double(),
                   Conductivity = col_double(),
                   Weather = col_character(),
                   Species = col_character(),
                   Mark = col_character(),
                   Catch = col_double(),
                   FL = col_double(),
                   Stage = col_character(),
                   Maturation = col_character(),
                   RaceByLength = col_character()
                 )) %>% 
  # Focusing on catch of LFS only
  mutate(Catch = ifelse(Species == "LFS", Catch, 0)) %>% 
  # Summing catch and averaging relevant environmental predictors to a daily
  # time step
  # Ignoring station and the number of tows in this demonstration
  group_by(Date) %>% 
  summarise(across(c(Volume, Secchi, DO, WaterTemp, Turbidity, Conductivity),
                   ~mean(.x, na.rm = T)),
            Catch = sum(Catch, na.rm = T))

# Defining the response variable: we can look at catch of longfin smelt here
# Since this is a multi-class classification example, catch will be transformed
# to a 3-level variable: no catch, low catch, high catch. High catch will be
# defined as any catch value greater than the median across the entire dataset
# Since catch is an ordered response here, the levels are specified
dataModel <- data %>% 
  mutate(Catch = factor(case_when(Catch == 0 ~ "noCatch",
                                  Catch > 0 & Catch < mean(Catch) ~ "lowCatch",
                                  Catch >= mean(Catch) ~ "highCatch"),
                        levels = c("noCatch", "lowCatch", "highCatch"))) %>% 
  # Removing date here as that will not be a predictor in the model
  select(-Date)

dim(dataModel)
# Number of missing datapoints:
colSums(is.na(dataModel))
# No missing catch data. Missing data in the environmental data but will leave 
# those as is in this demonstration

dataModel %>% 
  group_by(Catch) %>% 
  tally()
# Skewed distribution of our 3 classes: will use Kappa to account for this

# Quickly explore the dataset ---------------------------------------------

# This function is derived from the help page of ?pairs to show numberical cor
# Note that cor() here defaults to pearson and that we am asking for only
# complete pairs
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

# Correlation matrix
dataModel %>% 
  mutate(Catch = factor(Catch)) %>% 
  pairs(upper.panel = panel.cor, cex.labels = 2)
# Water temperature appears to be most correlated to catch of LFS, followed by
# Secchi and turbidity (of which are correlated to each other)

# Preprocessing the data --------------------------------------------------

# Use an 75/25 data split to train/test the model. This split is arbitrary.

set.seed(135)
splitIndex <- createDataPartition(y = dataModel$Catch, p = 0.75, list = F)

dataTrain <- dataModel[splitIndex, ]
dataTest <- dataModel[-splitIndex, ]

# Model training ----------------------------------------------------------

# Defining the seeds so that each train via CV can be reproduced
# The length of the list is the number of total reamples you want to use + 
# 1 final resample for the "optimal" tune across all your resamples
# Here, we will use a repeated CV approach with 10 folds and 3 repeats
seedList <- vector("list", 31)
set.seed(135)
for (i in 1:31) seedList[[i]] <- sample.int(n = 1000000, 1)

# First, fit the model against every feature (predictor) in the dataset
# The evaluation metric will be "Kappa" here as the response is not balanced
modelFullGini <- train(Catch ~ .,
                       data = dataTrain,
                       method = "rpart",
                       na.action = na.rpart,
                       trControl = trainControl(method = "repeatedcv",
                                                repeats = 3, 
                                                seeds = seedList,
                                                savePredictions = "final",
                                                classProb = T),
                       tuneLength = 10,
                       # # Can specify the specific tuning values if you do not
                       # # want to use tuneLength. This involves using tuneGrid
                       # # and the expand.grid() function for all tunable
                       # # parameters of your model, each specified as a column in
                       # # the grid. `rpart` only has one tunable
                       # # parameter, `cp`. The values below are what 
                       # # tuneLength = 10 picked for those that want to replicate it
                       # tuneGrid = expand.grid(
                       #   cp = c(0, 0.006001372, 0.012002743, 0.018004115,
                       #          0.024005487, 0.030006859, 0.036008230, 
                       #          0.042009602, 0.048010974, 0.054012346)
                       #   ),
                       metric = "Kappa")

# This variant uses entropy to make node splits
modelFullInfGain <- train(Catch ~ .,
                          data = dataTrain,
                          method = "rpart",
                          na.action = na.rpart,
                          parms = list(split = "information"),
                          trControl = trainControl(method = "repeatedcv",
                                                   repeats = 3, 
                                                   seeds = seedList,
                                                   savePredictions = "final",
                                                   classProb = T),
                          tuneLength = 10,
                          metric = "Kappa")

# For this dataset and application, the Gini selection measure performed
# slightly better in terms of Kappa

# Visualizing the tuning process for cp:
# In caret, the tuning process can be visualized with plot or ggplot
plot(modelFullGini)
# ggplot(modelFullGini) # (base layer only)
plot(modelFullInfGain)
# In this version of rpart in caret, the only tunable parameter is the cp
# For an exhaustive list of what parameters are available for tuning of a model
# in caret: https://topepo.github.io/caret/train-models-by-tag.html

# We see from the plot that our tune did find a maximum along its path and is
# appropriate

# Model performance evaluation --------------------------------------------

# The confusionMatrix function from caret is a useful tool to quickly understand
# the performance of a model. The function provides a confusion matrix and
# various evaluation statistics of the model predictions

# Across the training dataset:
# Gini model:
confusionMatrix(modelFullGini)
# InfGain model:
confusionMatrix(modelFullInfGain)
# Interestingly, the model did not predict any lowCatch class at all, in either
# variant.

# Across the testing dataset:
# Gini model:
confusionMatrix(predict(modelFullGini, dataTest, na.action = na.rpart), 
                dataTest$Catch)
# InfGain model:
confusionMatrix(predict(modelFullInfGain, dataTest, na.action = na.rpart), 
                dataTest$Catch)
# Both models classified similarly for setosa and versicolor. For virginia

# Compare model performances across the two selection measures variants in the
# training and testing datasets
data.frame(
  variant = c("Gini", "Gini", "infGain", "infGain"),
  dataSet = c("Train", "Test", "Train", "Test"),
  cpFinal = c(rep(filter(modelFullGini$results, 
                         Kappa == max(Kappa))[["cp"]], 2),
              rep(filter(modelFullInfGain$results, 
                         Kappa == max(Kappa))[["cp"]], 2)),
  Kappa = c(filter(modelFullGini$results, 
                   Kappa == max(Kappa))[["Kappa"]],
            confusionMatrix(predict(modelFullGini, newdata = dataTest,
                                    na.action = na.rpart), 
                            dataTest$Catch)$overall[["Kappa"]],
            # InfGain
            filter(modelFullInfGain$results, 
                   Kappa == max(Kappa))[["Kappa"]],
            confusionMatrix(predict(modelFullInfGain, newdata = dataTest,
                                    na.action = na.rpart), 
                            dataTest$Catch)$overall[["Kappa"]])
) %>%
  ggplot(aes(x = as.numeric(factor(dataSet, levels = c("Train", "Test"))),
             y = Kappa, 
             shape = variant,
             color = variant)) +
  geom_point(size = 4) +
  geom_line(size = 1.1) +
  scale_x_continuous(limits = c(0.75, 2.25),
                     breaks = c(1, 2),
                     labels = c("Train", "Test")) +
  labs(x = "Dataset",
       color = "Variant",
       shape = "Variant") +
  theme_bw(base_size = 24)
# The difference in Kappa between the training and testing datasets for the
# infGain variant is less than in the Gini variant. A smaller difference
# in the evaluation metrics between the training and testing datasets is 
# preferred; this indicates greater performance stability of the trained model 
# on unseen data.

# Variable importance -----------------------------------------------------
# Ranking order between both variants are the same
varImp(modelFullGini)
varImp(modelFullInfGain)

# Can get information on the relative error of the model as cp decreased:
# Other information as well such as which variables were used to construct
printcp(modelFullGini$finalModel)
printcp(modelFullInfGain$finalModel)

# Information on the node splits and the make up of the surrogate trees at each
# decision node
summary(modelFullGini$finalModel)
summary(modelFullInfGain$finalModel)

# Relevant visualizations -------------------------------------------------
# Visualizing the tree:
# Trees visualizations default to the T state splitting to the left of the tree
# Each leaf (terminal node) is labeled with a predicted class
# Values within each leaf depicts the distribution of datapoints captured in
# that leaf according to their true class (noCatch/lowCatch/highCatch)
par(mfrow = c(1, 2), mar = rep(0.1, 4))
# Gini version
plot(modelFullGini$finalModel, margin = 0.05)
text(modelFullGini$finalModel, use.n = T)
# IG version
plot(modelFullInfGain$finalModel, margin = 0.05)
text(modelFullInfGain$finalModel, use.n = T)
# Although the variants produced models with the same performance, you can see
# from the make up of the tree that classification of versicolor and virginica
# are slightly different between the two models. This slight difference should
# be considered by the modeler on which selection measure to use.

# Another way to visualize the tree:
# Color is added to this visualization, with intensity directly proportional 
# to the purity of the node (# predicted vs # observed)
rpart.plot(modelFullGini$finalModel, extra = 1)
rpart.plot(modelFullInfGain$finalModel, extra = 1)

# Visualizing the variable importance

print(plot(varImp(modelFullGini)), split = c(1, 1, 2, 1), more = T) 
print(plot(varImp(modelFullInfGain)), split = c(2, 1, 2, 1), more = F)

# Conclusion --------------------------------------------------------------
# We attempted to use a classification tree to describe the number of catch of
# LFS in the Chipps Island dataset. Catch was classified into three classes for
# this demonstration: no catch, low catch (catch < overall mean catch), and
# high catch (catch >= overall mean catch). Data was split into a training and
# a testing dataset via a 75/25 split. Two selection measures were tested
# to grow the tree, Gini and information gain. A 10-fold 3 repeat cross 
# validation approach was used to tune the cp metric. Model performance in the
# training and testing datasets across both selection measures were compared
# and the information gain variant preferred.

# Overall, high catch of LFS at Chipps is more likely when water temperature is
# low and water clarity is low. Curiously, the model did not predict much or 
# any low-catch in the training and testing dataset. This obviously 
# needs to be addressed. There are various ways to potentially improve model 
# fit: 1) better define the response variable as the distinction between the 
# three classes currently may not be precise enough to find relationships 
# between the predictor and response; 2) remove highly correlated predictors 
# (keep secchi and remove turbidity); 3) explore imputation for missing
# data in each predictor if appropriate; 4) add additional predictors (e.g., 
# day of year, although correlation to water temperature will need to be 
# explored); and 4) explore additional decision tree models (e.g., bagging, 
# boosting, ensemble).

# Misc code for markdown file ---------------------------------------------
# Quick model of the IRIS dataset; the entire dataset is used here to train
modelFullGiniIRIS <- train(Species ~ .,
                           data = iris,
                           trControl = trainControl(method = "repeatedcv",
                                                    repeats = 3, 
                                                    p = 0.75,
                                                    seeds = seedList,
                                                    verboseIter = F,
                                                    savePredictions = "final"),
                           tuneLength = 10,
                           method = "rpart",
                           na.action = na.rpart,
                           metric = "Kappa")

par(mfrow = c(1, 2), mar = rep(5, 4, 4, 2))
# Visualization of the decision tree
rpart.plot(modelFullGiniIRIS$finalModel)

# A decision tree partitions the feature space into sets of rectangles
plot(iris$Petal.Length, iris$Petal.Width, 
     xlab = "Petal Length",
     ylab = "Petal Width",
     pch = 20,
     cex = 2,
     col=ifelse(iris$Species == "setosa", "#FB6A4A", 
                ifelse(iris$Species == "virginica", "#74C476", 
                       ifelse(iris$Species ==  "versicolor",
                              "#999999", "white"))))
abline(v = 2.5)
abline(h = 1.8)
text(x = 1.6, y = 1.5, "setosa", col = "#FB6A4A", cex = 2)
text(x = 5, y = 0.5, "versicolor", col = "#999999", cex = 2)
text(x = 3.8, y = 2.3, "virginica", col = "#74C476", cex = 2)

# Saving RData file
save.image("ClassificationTrees.RData")
