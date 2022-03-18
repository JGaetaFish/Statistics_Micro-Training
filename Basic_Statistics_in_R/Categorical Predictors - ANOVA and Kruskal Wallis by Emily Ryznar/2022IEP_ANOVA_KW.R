#IEP Workshop Training: ANOVA, Kruskal-Wallis, and Multiple Comparisons

#By: Emily Ryznar

#3/23/2022


##QUESTION: Are Chinook fork length similar among life history stages in 2004?

#############################################################################################
#Import and and prepare data                                                                #
#############################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Import and view data
data<- read.csv(file="LI_seine.csv", header=T) ##import data

head(data) #view data

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Pull out sampling years into new column
library(lubridate) #load package

data$SampleDate<-mdy(data$SampleDate) #switching date to standard format for lubridate
data$SampleYear<-year(data$SampleDate) #pulling out sample year from sample date, store in new column
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Subset data for only Chinook salmon stages 2:5 in 2004 and label stages

ch_2004<-subset(data,data$StageCode!="1" & data$CommonName=="Chinook Salmon"
                & data$SampleYear=="2004") #subset data

labels = data.frame(StageCode = 2:5,
                    stage = c("Fry", "Parr", "Silvery parr", "Smolt")) #generate stage labels

ch_2004<-merge(ch_2004, labels, all.x=TRUE) #merge labels data frame with salmon data frame
head(ch_2004)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Plot raw data by stage, visually inspect for outliers (though can't remove because no justification)
dev.off()
plot(ForkLength ~ jitter(StageCode, factor=1/3), pch=20,
     col=gray(0.25,0.5), data=ch_2004, xaxt="n", las=1,
     ylab="Fork Length (mm)", xlab="Life Stage")
axis(side = 1, at=2:5,
     labels = c("Fry", "Parr", "Silvery Parr", "Smolt")) #plot data

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Check for missing values

unique(is.na(ch_2004$ForkLength)) #false=none!

#############################################################################################
# 1-Factor ANOVA                                                                            #
#############################################################################################

##QUESTION: Are Chinook fork lengths similar among life history stages in 2004?

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Make sure "stage" is a factor; run 1-Factor ANOVA on untransformed data and untransformed data

ch_2004$stage <- as.factor(ch_2004$stage) #stage as factor

un_test<-lm(ForkLength~stage, ch_2004) #ANOVA as a linear model with untransformed data

test <- lm((1/ForkLength) ~ stage, ch_2004) #ANOVA as a linear model with 1/x transformed data

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Visually inspect fitted model residuals (using untransformed and transformed data)
##for normality and equal variances (ANOVA assumptions) 

par(mfrow = c(2, 2), mar = c(4.1, 4.1, 0.5, 0.5)) #specify plotting setup

hist(un_test$residuals, main = NULL, cex.axis = 0.9, cex.lab = 0.9, ylab = "", xlab = "")
title(ylab = "Frequency", xlab = "un_test$residuals", line = 2, cex.lab = 0.9)
box(which = "plot")
text(-24, 90, "(a)") ##histogram for untransformed normality

plot(un_test$residuals ~ ch_2004$stage, cex.axis = 0.9, cex.lab = 0.9, ylab = "", xlab = "")
title(ylab = "un_test$residuals", xlab = "Stage", line = 2, cex.lab = 0.9)
text(0.6, 18, "(b)") ##boxplot for untransformed variance

hist(test$residuals, main = NULL, cex.axis = 0.9, cex.lab = 0.9, ylab = "", xlab = "")
title(ylab = "Frequency", xlab = "test$residuals", line = 2, cex.lab = 0.9)
box(which = "plot")
text(-0.0075, 103, "(c)") #histogram for 1/x transformed normality

plot(test$residuals ~ ch_2004$stage, cex.axis = 0.9, cex.lab = 0.9, ylab = "", xlab = "")
title(ylab = "test$residuals", xlab = "Stage", line = 2, cex.lab = 0.9)
text(0.6, 0.0073, "(d)") #boxplot for 1/x transformed variance

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Statistically assess fitted model residuals (using untransformed and transformed data)
##for normality and equal variances (ANOVA assumptions)

shapiro.test(un_test$residuals) #untransformed normality test

shapiro.test(test$residuals) #transformed normality test

library(car) #load "car" package for leveneTest() function below

bartlett.test(un_test$residuals ~ stage, ch_2004) #untransformed variance test
leveneTest(un_test$residuals ~ stage, ch_2004) #untransformed variance test

bartlett.test(test$residuals ~ stage, ch_2004) #transformed variance test
leveneTest(test$residuals ~ stage, ch_2004) #transformed variance test

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##View ANOVA model summary using 1/x transformed data 

summary(test) #view summary of ANOVA output

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Run multiple comparisons
library(multcomp) #load 'multcomp' package for multiple comparisons test

tuk <- glht(test, linfct = mcp(stage = "Tukey")) #assess pairwise differences using Tukey HSD

summary(tuk) #view pairwise summary

cld(tuk) #view compact letter display of pairwise summary

#############################################################################################
# Kruskal-Wallis                                                                            #
#############################################################################################

##QUESTION: Are Chinook fork length similar among life history stages in 2004?

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Visually inspect for similar fork length distribution among stages (KW assumption)

par(mfrow=c(2,2), mar=c(4.5,4.5,0.5,0.5)) #specify plotting window

hist(ch_2004$ForkLength[which(ch_2004$stage=="Fry")], main=NULL, cex.axis=0.9, 
     cex.lab=0.9, ylab=NULL, xlab=NULL)
title(ylab="Frequency", xlab="Fry fork length (mm)", line=2, cex.lab=0.9)
box(which = "plot")
text(34.5, 25, "(a)") #histogram for fry fork length distribution

hist(ch_2004$ForkLength[which(ch_2004$stage=="Parr")], main=NULL, cex.axis=0.9, 
     cex.lab=0.9, ylab=NULL, xlab=NULL)
title(ylab="Frequency", xlab="Parr fork length (mm)", line=2, cex.lab=0.9)
box(which = "plot")
text(36, 21, "(b)") #histogram for parr fork length distribution

hist(ch_2004$ForkLength[which(ch_2004$stage=="Silvery parr")], main=NULL, cex.main=0.9, 
     cex.axis=0.9, cex.lab=0.9, ylab=NULL, xlab=NULL)
title(ylab="Frequency", xlab="Silvery parr fork length (mm)", line=2, cex.lab=0.9)
box(which = "plot")
text(41, 23, "(c)") #histogram for silvery parr length distribution

hist(ch_2004$ForkLength[which(ch_2004$stage=="Smolt")], main=NULL, cex.axis=0.9, 
     cex.lab=0.9, ylab=NULL, xlab=NULL)
title(ylab="Frequency", xlab="Smolt fork length (mm)", line=2, cex.lab=0.9)
box(which = "plot")
text(56, 3.7, "(d)") #histogram for smolt fork length distribution

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Run Kruskal-Wallis analysis and view results

library(FSA) #load "FSA" package for kruskal.test()

kw_test <- kruskal.test(ForkLength ~ stage, ch_2004) #run Kruskal-Wallis test

kw_test #view results

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Run pairwise Kruskal-Wallis comparison

library(stats) #load "stats" package for pairwise.wilcox.test()

pairwise.wilcox.test(ch_2004$ForkLength, ch_2004$stage, 
                     p.adjust.method = "BH") #run pairwise Wilcox test with BH correction

#############################################################################################
# Model presentation                                                                        #
#############################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Generate boxplot of fork length by stage with means plotted as blue dots and Tukey post-hoc labels added
dev.off()
boxplot(ch_2004$ForkLength~ch_2004$stage, ylim=c(30, 95), ylab="Fork length (mm)", xlab="Stage")
points(1, mean(ch_2004$ForkLength[which(ch_2004$stage=="Fry")]), col="blue", pch=19, cex=1.25) #add mean
points(2, mean(ch_2004$ForkLength[which(ch_2004$stage=="Parr")]), col="blue", pch=19, cex=1.25)
points(3,mean(ch_2004$ForkLength[which(ch_2004$stage=="Silvery parr")]), col="blue", pch=19, cex=1.25)
points(4,mean(ch_2004$ForkLength[which(ch_2004$stage=="Smolt")]), col="blue", pch=19, cex=1.25)
text(1, max(ch_2004$ForkLength[which(ch_2004$stage=="Fry")])+3.5, "a") #add post hoc labels
text(2, max(ch_2004$ForkLength[which(ch_2004$stage=="Parr")])+3.5, "b")
text(3, max(ch_2004$ForkLength[which(ch_2004$stage=="Silvery parr")])+3.5, "c")
text(4, max(ch_2004$ForkLength[which(ch_2004$stage=="Smolt")])+3.5, "d")
















