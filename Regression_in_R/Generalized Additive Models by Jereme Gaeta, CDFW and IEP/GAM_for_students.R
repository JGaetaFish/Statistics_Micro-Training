##############################################################
# Generalized Additive Modeling
# IEP Micro-Training
# Created by Jereme W. Gaeta, IEP & CDFW
# email: Jereme.Gaeta@Wildlife.ca.gov
##############################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Reset your environment:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# rm(list=ls())
# graphics.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Load Libraries and custom function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

detach(package:gam)
library(mgcv)
library(itsadug)
library(lubridate)
library(data.table)

#~ A custom function to calculate day of water year
water.day = function(x, start.month = 10L){
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1L)
  as.integer(x - start.date + 1L)
}

#~ A custom function to add a label to the same location on every panel
plot_label = function(lab="(a)", x_prop=0.08, y_prop=0.92,
                      font_type=2, fcex=1.15, usr=par('usr')){
  x_val = usr[1]+(usr[2]-usr[1])*x_prop
  y_val = usr[3]+(usr[4]-usr[3])*y_prop
  text(x = x_val, y = y_val, labels = lab, font = font_type, cex=fcex)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Load and explore the data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df = read.csv(file = "Yolo_tow_drain_drift_cpue.csv", header=T)

head(df)


at_dates = water.day(as.Date(c("2021-10-1", "2021-12-1", "2022-2-1", "2022-4-1", 
                               "2022-6-1", "2022-8-1", "2022-9-30"), format="%Y-%m-%d"))
lab_dates = format(as.Date(c("2021-10-1", "2021-12-1", "2022-2-1", "2022-4-1", 
                             "2022-6-1", "2022-8-1", "2022-9-30"),
                           format="%Y-%m-%d"), format="%h-%d")

par(mfrow=c(1,2), mar=c(4,4.5,4,1)+0.1, oma=c(0,0,0,0))
hist(df$cpue, main="", ylab="", xlab="", las=1)
mtext(text = "CPUE", side = 1, line = 2.25)
mtext(text = "Frequency", side = 2, line = 3)
box(which="plot")
plot_label(lab = "(a)", x_prop = 0.92)

hist(log(df$cpue), main="", ylab="", xlab="", las=1)
mtext(text = expression(log[e]~"(CPUE)"), side = 1, line = 2.25)
mtext(text = "Frequency", side = 2, line = 3)
box(which="plot")
plot_label(lab = "(b)", x_prop = 0.92)



at_dates = water.day(as.Date(c("2021-10-1", "2021-12-1", "2022-2-1", "2022-4-1", 
                               "2022-6-1", "2022-8-1"), format="%Y-%m-%d"))
lab_dates = format(as.Date(c("2021-10-1", "2021-12-1", "2022-2-1", "2022-4-1", 
                             "2022-6-1", "2022-8-1"), format="%Y-%m-%d"),
                   format="%h-%d")
lm_mod = lm(log(cpue) ~ water_doy, data = df)
lm_dat = data.frame(water_doy = seq(from=min(df$water_doy),
                                    to=max(df$water_doy),
                                    length.out=100))
lm_pred = predict(object = lm_mod, newdata = lm_dat)
resid_supsmu=supsmu(x=df$water_doy, y=residuals(lm_mod))

par(mfrow=c(1,3), mar=c(4,4.5,4,1)+0.1, oma=c(0,0,0,0))
plot(df$cpue ~ df$water_doy, pch=20, col=gray(0.1,0.2),
     ylab="", xlab="", las=1)
axis(side = 3, at = at_dates, labels = lab_dates, las=2, cex=0.8)
mtext(text = "CPUE", side = 2, line = 3.25)
mtext(text = "Day of Water Year", side = 1, line = 2.25)
plot_label(lab = "(a)")

plot(log(df$cpue) ~ df$water_doy, pch=20, col=gray(0.1,0.2),
     ylab="", xlab="", las=1)
axis(side = 3, at = at_dates, labels = lab_dates, las=2, cex=0.8)
mtext(text = expression(log[e]~"(CPUE)"), side = 2, line = 2.25)
mtext(text = "Day of Water Year", side = 1, line = 2.5)
plot_label(lab = "(b)")
lines(lm_pred ~ lm_dat$water_doy, lwd=2)

plot(residuals(lm_mod) ~ df$water_doy, pch=20, col=gray(0.1,0.2),
     ylab="", xlab="", las=1)
axis(side = 3, at = at_dates, labels = lab_dates, las=2, cex=0.8)
mtext(text = "Model Residual", side = 2, line = 2.25)
mtext(text = "Day of Water Year", side = 1, line = 2.5)
abline(h=0, col=gray(0.5,0.5), lwd=3)
lines(x = resid_supsmu$x, y = resid_supsmu$y, col = "blue")
plot_label(lab = "(c)")




##############################################################
# Single Smoother Generalized Additive Model (GAM)
##############################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Build the GAM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#The plots above indicate we should log-e transform the data
df$ln_cpue = log(df$cpue)

mod = gam(ln_cpue ~ s(water_doy, bs = "cr"), data = df)
summary(mod)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Default GAM plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

par(mar=c(4,4,0.25,0.25)+0.1, oma=c(0,0,0,0))
plot(mod)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Evaluate model fit:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

layout(mat = matrix(c(1,2,3,3), 2, 2, byrow=TRUE))
par(mar=c(4.5,4.5,1.5,1.5)+0.1, oma=rep(0,4))
plot(df$ln_cpue ~ predict(mod), las=1, pch=20, col=gray(0.2,0.2))
abline(0,1)
plot_label(lab = "(a)")
hist(residuals(mod), las=1, main=""); box(which="plot")
plot_label(lab = "(b)")
plot(residuals(mod) ~ df$water_doy, las=1, pch=20, col=gray(0.2,0.2))
abline(h=0)
plot_label(lab = "(c)")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Visualize model uncertainty:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Determine the confidence interval to present in your figure and determine the
# t-distribution value of which to multiply the standard error
interval_value = 0.95
t_distribution_probability = 1.0 - (1-interval_value)/2
CI_factor = qt(t_distribution_probability, Inf)

#~ Generate a data frame with the sequence from the minimum to the maximum
# observed value of the predictor variable: day of water year
mod_doy_range = range(na.omit(df$water_doy))
new_dat_1 = data.frame(water_doy = seq(from = mod_doy_range[1],
                                       to = mod_doy_range[2], by=1))

#~ Use the *predict()* function to estimate model fit and standard error
preds = predict(mod, newdata = new_dat_1, se.fit = TRUE)

# Multiply the standard error fit by the t-distribution value to estimate 
# the upper and lower confidence interval limits
fit = preds$fit
upper = fit-CI_factor*preds$se.fit
lower = fit+CI_factor*preds$se.fit


#~~   Visualize the model

# Plot the raw data and make the axes and labels aesthetically pleasing
par(mfrow=c(1,1),mar=c(2.5, 4.5, 2,4.5)+0.1, oma=c(1.5,0,0.5,0))
plot(df$ln_cpue ~ df$water_doy, pch=20, col=gray(0.1,0.2),
     las=1, ylab="", xlab="")
mtext(text = "CPUE", side = 4, line = 3)
mtext(text = expression(log[e]~"(CPUE)"), side = 2, line = 2.5)
y_labs = c(0.0001,0.0003,0.001, 0.003,0.01,0.033,
           0.1,0.33, 1, 3, 8, 20, 55, 150)[c(T,F)]
axis(side = 4, at = log(y_labs), labels = y_labs, las=1)
mtext(text = "Day of Water Year", side = 1, line = 2.25)

# Use the polygon function to add the CI
polygon(x = c(new_dat_1$water_doy, max(new_dat_1$water_doy), 
              rev(new_dat_1$water_doy), new_dat_1$water_doy[1]),
        y = c(lower, upper[length(upper)],
              rev(upper), lower[1]), border=NA,
        col=rgb(20,200,20,alpha=150,maxColorValue=255))

# Use the lines() function to add the mean model prediction
lines(x = new_dat_1$water_doy, y = fit, lwd=2,
      col=rgb(20,150,20,alpha=255,maxColorValue=255))


##############################################################
# Multiple Smoother Generalized Additive Model (GAM)
##############################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Explore the data: Yolo Bypass Inundation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

at_dates = water.day(as.Date(c("2021-10-1", "2021-12-1", "2022-2-1", "2022-4-1", 
                               "2022-6-1", "2022-8-1", "2022-10-1"),
                             format="%Y-%m-%d"))
lab_dates = format(as.Date(c("2021-10-1", "2021-12-1", "2022-2-1", "2022-4-1", 
                             "2022-6-1", "2022-8-1", "2022-10-1"),
                           format="%Y-%m-%d"), format="%h-%d")
poly_col = c("firebrick3", "dodgerblue4")
col_mat = as.matrix(col2rgb(poly_col))

layout(mat = matrix(c(1,2,2), 1, 3, byrow=TRUE))
par(oma=c(0,0,0,2))
boxplot(df$ln_cpue ~ df$fInun, ylab="",
        names = NA, whisklty=1,
        col=rgb(red = col_mat[1,],green = col_mat[2,],
                blue = col_mat[3,], alpha=200,
                maxColorValue = 255), las=1,
        xlab="")
axis(side = 1, at=c(1,2), labels = c("", "Inundated"), line = 0.25,
     tick=FALSE)
axis(side = 1, at=c(1,2), labels = c("Not\nInundated", ""),
     line = 0.75, tick = FALSE)
mtext(text = expression(log[e]~"(CPUE)"), side = 2, line = 2.5)
plot_label(lab = "(a)")

plot(df$ln_cpue ~ df$water_doy, type="n",
     las=1, ylab="", xlab="")
mtext(text = "CPUE", side = 4, line = 3)
mtext(text = expression(log[e]~"(CPUE)"), side = 2, line = 2.5)
y_labs = c(0.0001,0.0003,0.001, 0.003,0.01,0.033,
           0.1,0.33, 1, 3, 8, 20, 55, 150)[c(T,F)]
axis(side = 4, at = log(y_labs), labels = y_labs, las=1)
mtext(text = "Day of Water Year", side = 1, line = 2.25)
axis(side = 3, at = at_dates, labels = lab_dates)
plot_inun = c(0, 1)
for(i in 1:length(plot_inun)){
  sub = subset(df, df$fInun==plot_inun[i])
  points(sub$ln_cpue ~ sub$water_doy, pch=16,
         col=rgb(red = col_mat[1,i],green = col_mat[2,i],
                 blue = col_mat[3,i], alpha=200,
                 maxColorValue = 255))
}
plot_label(lab = "(b)")
legend("topright", legend = c("Inundated", "Not Inundated"), pch=16, 
       col=rev(rgb(red = col_mat[1,],green = col_mat[2,],
                   blue = col_mat[3,], alpha=200,
                   maxColorValue = 255)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Build the GAM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# make sure fInun is, in fact, a factor
df$fInun = as.factor(df$fInun)

#build the model
mod2 = gam(ln_cpue ~ fInun + s(water_doy, bs = "cr", by = fInun), 
           data = df)

summary(mod)
summary(mod2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Evaluate model fit:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inun_index_list = list(inun0 = which(df$fInun=="0"),
                       inun1 = which(df$fInun=="1"))

layout(mat = matrix(c(1,1,1,2,2,2,3,3,3,3,4,4), 2, 6, byrow=TRUE))
par(mar=c(4.5,4.5,1.5,1.5)+0.1, oma=rep(0,4))
plot(df$ln_cpue ~ predict(mod2), las=1, type="n")
for(i in 1:2){
  points(df$ln_cpue[inun_index_list[[i]]] ~
           predict(mod2)[inun_index_list[[i]]], pch=16,
         col=rgb(red = col_mat[1,i],green = col_mat[2,i],
                 blue = col_mat[3,i], alpha=200,
                 maxColorValue = 255))
}
abline(0,1)
plot_label(lab = "(a)")

hist(residuals(mod2), main="",
     col=rgb(red = col_mat[1,1],green = col_mat[2,1],
             blue = col_mat[3,1], alpha=255,
             maxColorValue = 255), 
     breaks = c(seq(from=-8, to=8, by=0.5)))
hist(residuals(mod2)[inun_index_list[[2]]], add = T,
     col=rgb(red = col_mat[1,2],green = col_mat[2,2],
             blue = col_mat[3,2], alpha=255,
             maxColorValue = 255), 
     breaks = c(seq(from=-8, to=8, by=0.5)))
box(which="plot")
plot_label(lab = "(b)")

plot(residuals(mod2) ~ df$water_doy, las=1, type="n")
for(i in 1:2){
  points(residuals(mod2)[inun_index_list[[i]]] ~
           df$water_doy[inun_index_list[[i]]], pch=16,
         col=rgb(red = col_mat[1,i],green = col_mat[2,i],
                 blue = col_mat[3,i], alpha=200,
                 maxColorValue = 255))
}
abline(h=0)
plot_label(lab = "(c)")

plot(residuals(mod2) ~ df$fInun, las=1, pch=20, 
     xlab="", whisklty=1,
     col=rgb(red = col_mat[1,],green = col_mat[2,],
             blue = col_mat[3,], alpha=255,
             maxColorValue = 255))
abline(h=0)
plot_label(lab = "(d)")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Visualize model uncertainty:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

prob_values = c(0.00001, 0.999999)
prob_mat = matrix(data = NA, nrow = 2, ncol = 2, byrow = TRUE,
                  dimnames = list(inundation = c(0,1),
                                  water_doy = paste0("Prob_", prob_values)))
for(i in 1:length(levels(df$fInun))){
  sub = subset(df, df$fInun==levels(df$fInun)[i])
  prob_mat[i,] = round(quantile(x = na.omit(sub$water_doy),
                                probs=prob_values),
                       digits=0)
}

inun_df = data.frame(water_doy = prob_mat[2,1]:prob_mat[2,2],
                     fInun = rep("1", times=length(prob_mat[2,1]:prob_mat[2,2])))
dry_df = data.frame(water_doy = prob_mat[1,1]:prob_mat[1,2],
                    fInun = rep("0", times=length(prob_mat[1,1]:prob_mat[1,2])))

inun_pred = predict(mod2, newdata = inun_df, se.fit = TRUE)
dry_pred = predict(mod2, newdata = dry_df, se.fit = TRUE)

inun_fit = inun_pred$fit
inun_upper = inun_fit-CI_factor*inun_pred$se.fit
inun_lower = inun_fit+CI_factor*inun_pred$se.fit

dry_fit = dry_pred$fit
dry_upper = dry_fit-CI_factor*dry_pred$se.fit
dry_lower = dry_fit+CI_factor*dry_pred$se.fit

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~   Visualize the model: multiple smoothers

par(mfrow=c(1,1),mar=c(2.5, 4.5, 2,4.5)+0.1, oma=c(1.5,0,0.5,0))
plot(df$ln_cpue ~ df$water_doy, type="n",
     las=1, ylab="", xlab="")
plot_inun = c(0, 1)
for(i in 1:length(plot_inun)){
  sub = subset(df, df$fInun==plot_inun[i])
  points(sub$ln_cpue ~ sub$water_doy, pch=16,
         col=rgb(red = col_mat[1,i],green = col_mat[2,i],
                 blue = col_mat[3,i], alpha=180,
                 maxColorValue = 255))
}
mtext(text = "CPUE", side = 4, line = 3)
mtext(text = expression(log[e]~"(CPUE)"), side = 2, line = 2.5)
y_labs = c(0.0001,0.0003,0.001, 0.003,0.01,0.033,
           0.1,0.33, 1, 3, 8, 20, 55, 150)[c(T,F)]
axis(side = 4, at = log(y_labs), labels = y_labs, las=1)
mtext(text = "Day of Water Year", side = 1, line = 2.25)
polygon(x = c(inun_df$water_doy, max(inun_df$water_doy),
              rev(inun_df$water_doy), inun_df$water_doy[1]),
        y = c(inun_lower, inun_upper[length(inun_upper)], 
              rev(inun_upper), inun_lower[1]), 
        col=rgb(red = col_mat[1,2],green = col_mat[2,2],
                blue = col_mat[3,2], alpha=200,
                maxColorValue = 255),
        border=NA)
lines(inun_fit ~ inun_df$water_doy, lwd=3,
      col="darkblue")
polygon(x = c(dry_df$water_doy, max(dry_df$water_doy),
              rev(dry_df$water_doy), dry_df$water_doy[1]),
        y = c(dry_lower, dry_upper[length(dry_upper)], 
              rev(dry_upper), dry_lower[1]), 
        col=rgb(red = col_mat[1,1],green = col_mat[2,1],
                blue = col_mat[3,1], alpha=180,
                maxColorValue = 255),
        border=NA)
lines(dry_fit ~ dry_df$water_doy, lwd=3,
      col="darkred")
axis(side = 3, at = at_dates, labels = lab_dates)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Visualize model uncertainty without data:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

prob_values = c(0.1, 0.9)
prob_mat = matrix(data = NA, nrow = 2, ncol = 2, byrow = TRUE,
                  dimnames = list(inundation = c(0,1),
                                  water_doy = paste0("Prob_", prob_values)))
for(i in 1:length(levels(df$fInun))){
  sub = subset(df, df$fInun==levels(df$fInun)[i])
  prob_mat[i,] = round(quantile(x = na.omit(sub$water_doy),
                                probs=prob_values),
                       digits=0)
}

inun_df = data.frame(water_doy = prob_mat[2,1]:prob_mat[2,2],
                     fInun = rep("1", times=
                                   length(prob_mat[2,1]:
                                            prob_mat[2,2])))
dry_df = data.frame(water_doy = prob_mat[1,1]:prob_mat[1,2],
                    fInun = rep("0", times=
                                  length(prob_mat[1,1]:
                                           prob_mat[1,2])))

inun_pred = predict(mod2, newdata = inun_df, se.fit = TRUE)
dry_pred = predict(mod2, newdata = dry_df, se.fit = TRUE)

inun_fit = inun_pred$fit
inun_upper = inun_fit-CI_factor*inun_pred$se.fit
inun_lower = inun_fit+CI_factor*inun_pred$se.fit

dry_fit = dry_pred$fit
dry_upper = dry_fit-CI_factor*dry_pred$se.fit
dry_lower = dry_fit+CI_factor*dry_pred$se.fit


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~   Visualize the model: multiple smoothers

par(mfrow=c(2,1),mar=c(2.5, 4.5, 2,4.5)+0.1, oma=c(1.5,0,0.5,0))
plot(df$ln_cpue ~ df$water_doy, type="n",
     las=1, ylab="", xlab="")
plot_inun = c(0, 1)
for(i in 1:length(plot_inun)){
  sub = subset(df, df$fInun==plot_inun[i])
  points(sub$ln_cpue ~ sub$water_doy, pch=16,
         col=rgb(red = col_mat[1,i],green = col_mat[2,i],
                 blue = col_mat[3,i], alpha=180,
                 maxColorValue = 255))
}
plot_label(lab = "(a)")
mtext(text = "CPUE", side = 4, line = 3)
mtext(text = expression(log[e]~"(CPUE)"), side = 2, line = 2.5)
y_labs = c(0.0001,0.0003,0.001, 0.003,0.01,0.033,
           0.1,0.33, 1, 3, 8, 20, 55, 150)[c(T,F)]
axis(side = 4, at = log(y_labs), labels = y_labs, las=1)
polygon(x = c(inun_df$water_doy, max(inun_df$water_doy),
              rev(inun_df$water_doy), inun_df$water_doy[1]),
        y = c(inun_lower, inun_upper[length(inun_upper)], 
              rev(inun_upper), inun_lower[1]), 
        col=rgb(red = col_mat[1,2],green = col_mat[2,2],
                blue = col_mat[3,2], alpha=200,
                maxColorValue = 255),
        border=NA)
lines(inun_fit ~ inun_df$water_doy, lwd=3,
      col="darkblue")
polygon(x = c(dry_df$water_doy, max(dry_df$water_doy),
              rev(dry_df$water_doy), dry_df$water_doy[1]),
        y = c(dry_lower, dry_upper[length(dry_upper)], 
              rev(dry_upper), dry_lower[1]), 
        col=rgb(red = col_mat[1,1],green = col_mat[2,1],
                blue = col_mat[3,1], alpha=180,
                maxColorValue = 255),
        border=NA)
lines(dry_fit ~ dry_df$water_doy, lwd=3,
      col="darkred")
axis(side = 3, at = at_dates, labels = lab_dates)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(exp(range(inun_upper, inun_lower, dry_upper,dry_lower)) ~
       range(dry_df$water_doy), type = "n", las=1,
     ylab="", xlab="", xlim = range(df$water_doy))
plot_label(lab = "(b)")
mtext(text = "Day of Water Year", side = 1, line = 2.5)
mtext(text = "CPUE", side = 2, line = 2.75)
polygon(x = c(inun_df$water_doy, max(inun_df$water_doy),
              rev(inun_df$water_doy), inun_df$water_doy[1]),
        y = exp(c(inun_lower, inun_upper[length(inun_upper)], 
                  rev(inun_upper), inun_lower[1])), 
        col=rgb(red = col_mat[1,2],green = col_mat[2,2],
                blue = col_mat[3,2], alpha=200,
                maxColorValue = 255),
        border=NA)
lines(exp(inun_fit) ~ inun_df$water_doy, lwd=3,
      col="darkblue")
polygon(x = c(dry_df$water_doy, max(dry_df$water_doy),
              rev(dry_df$water_doy), dry_df$water_doy[1]),
        y = exp(c(dry_lower, dry_upper[length(dry_upper)], 
                  rev(dry_upper), dry_lower[1])), 
        col=rgb(red = col_mat[1,1],green = col_mat[2,1],
                blue = col_mat[3,1], alpha=180,
                maxColorValue = 255),
        border=NA)
lines(exp(dry_fit) ~ dry_df$water_doy, lwd=3,
      col="darkred")
rug(df$water_doy[which(df$fInun=="0")], ticksize = 0.03,
    side = 1, lwd = 0.75,
    col=rgb(red = col_mat[1,1],green = col_mat[2,1],
            blue = col_mat[3,1], alpha=180,
            maxColorValue = 255))
rug(df$water_doy[which(df$fInun=="1")], ticksize = 0.03,
    side = 3, lwd = 0.75,
    col=rgb(red = col_mat[1,2],green = col_mat[2,2],
            blue = col_mat[3,2], alpha=200,
            maxColorValue = 255))
axis(side = 3, at = at_dates, labels = lab_dates)



##############################################################
# GAMs and Statistical Inference
##############################################################

# Determine the confidence interval to present in your figure and determine the
# t-distribution value of which to multiply the standard error
interval_value = 0.835
t_distribution_probability = 1.0 - (1-interval_value)/2
CI_factor = qt(t_distribution_probability, Inf)


prob_values = c(0.1, 0.9)
prob_mat = matrix(data = NA, nrow = 2, ncol = 2, byrow = TRUE,
                  dimnames = list(inundation = c(0,1),
                                  water_doy = paste0("Prob_", prob_values)))
for(i in 1:length(levels(df$fInun))){
  sub = subset(df, df$fInun==levels(df$fInun)[i])
  prob_mat[i,] = round(quantile(x = na.omit(sub$water_doy),
                                probs=prob_values),
                       digits=0)
}

inun_df = data.frame(water_doy = prob_mat[2,1]:prob_mat[2,2],
                     fInun = rep("1", times=
                                   length(prob_mat[2,1]:
                                            prob_mat[2,2])))
dry_df = data.frame(water_doy = prob_mat[1,1]:prob_mat[1,2],
                    fInun = rep("0", times=
                                  length(prob_mat[1,1]:
                                           prob_mat[1,2])))

inun_pred = predict(mod2, newdata = inun_df, se.fit = TRUE)
dry_pred = predict(mod2, newdata = dry_df, se.fit = TRUE)

inun_fit = inun_pred$fit
inun_upper = inun_fit-CI_factor*inun_pred$se.fit
inun_lower = inun_fit+CI_factor*inun_pred$se.fit

dry_fit = dry_pred$fit
dry_upper = dry_fit-CI_factor*dry_pred$se.fit
dry_lower = dry_fit+CI_factor*dry_pred$se.fit



#~~~~  Stat difference: bottom of the heap:
# https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/

smooth_diff <- function(model, newdata, f1, f2, var, alpha = 0.05,
                        unconditional = FALSE) {
  xp <- predict(model, newdata = newdata, type = 'lpmatrix')
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
  crit <- qt(alpha/2, df.residual(model), lower.tail = TRUE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  data.frame(pair = paste(f1, f2, sep = '-'),
             diff = dif,
             se = se,
             upper = upr,
             lower = lwr)
}
pdat <- expand.grid(water_doy = prob_mat[2,1]:prob_mat[2,2],
                    fInun=c("0", "1"))
comp1 <- smooth_diff(mod2, pdat, '0', '1', 'fInun')
comp <- cbind(water_doy = prob_mat[2,1]:prob_mat[2,2],
              comp1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

par(mfrow=c(2,1),mar=c(2.5, 4.5, 2,4.5)+0.1, oma=c(1.5,0,0.5,0))
plot(range(inun_upper, inun_lower, dry_upper,dry_lower) ~
       range(df$water_doy), type = "n", las=1,
     ylab="", xlab="")
plot_label(lab = "(a)")
mtext(text = "CPUE", side = 4, line = 3)
mtext(text = expression(log[e]~"(CPUE)"), side = 2, line = 2.75)
y_labs = c(0.01, 0.02, 0.03, 0.05, 0.08, 0.14, 0.22, 0.37)
axis(side = 4, at = log(y_labs), labels = y_labs, las=1)
polygon(x = c(inun_df$water_doy, max(inun_df$water_doy),
              rev(inun_df$water_doy), inun_df$water_doy[1]),
        y = c(inun_lower, inun_upper[length(inun_upper)], 
              rev(inun_upper), inun_lower[1]), 
        col=rgb(red = col_mat[1,2],green = col_mat[2,2],
                blue = col_mat[3,2], alpha=200,
                maxColorValue = 255),
        border=NA)
lines(inun_fit ~ inun_df$water_doy, lwd=3,
      col="darkblue")
polygon(x = c(dry_df$water_doy, max(dry_df$water_doy),
              rev(dry_df$water_doy), dry_df$water_doy[1]),
        y = c(dry_lower, dry_upper[length(dry_upper)], 
              rev(dry_upper), dry_lower[1]), 
        col=rgb(red = col_mat[1,1],green = col_mat[2,1],
                blue = col_mat[3,1], alpha=180,
                maxColorValue = 255),
        border=NA)
lines(dry_fit ~ dry_df$water_doy, lwd=3,
      col="darkred")
rug(df$water_doy[which(df$fInun=="0")], ticksize = 0.03,
    side = 1, lwd = 0.75,
    col=rgb(red = col_mat[1,1],green = col_mat[2,1],
            blue = col_mat[3,1], alpha=180,
            maxColorValue = 255))
rug(df$water_doy[which(df$fInun=="1")], ticksize = 0.03,
    side = 3, lwd = 0.75,
    col=rgb(red = col_mat[1,2],green = col_mat[2,2],
            blue = col_mat[3,2], alpha=200,
            maxColorValue = 255))
axis(side = 3, at = at_dates, labels = lab_dates)

plot(range(comp$upper, comp$lower) ~ range(df$water_doy),
     type = "n", las=1, ylab="", xlab="")
mtext(text = "Difference in CPUE Trend", side = 2, line = 2.75)
mtext(text = "Day of Water Year", side = 1, line = 2.5)
plot_label(lab = "(b)")
abline(h=0, lty=1, lwd=2)
polygon(x = c(comp$water_doy, max(comp$water_doy), rev(comp$water_doy),
              comp$water_doy[1]),
        y = c(comp$lower, comp$upper[dim(comp)[1]], rev(comp$upper),
              comp$lower[1]), col=gray(0.5,0.5), border=NA)
lines(comp$diff ~ comp$water_doy, lwd=3, col = gray(0.5,1))

##############################################################
# Autocorrelation Assessment
##############################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Data autocorrelation assessment:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# make sure the data are sorted temporally before we begin
df=df[order(df$water_year, df$water_doy, decreasing = FALSE),]

## lines beyond the blue dashed line are lags 
par(mfrow=c(2,1), mar=c(3,4.5,1,1)+0.1, oma=c(1,0,0,0))
acf(df$ln_cpue, main="", las=1)
plot_label(lab = "(a)")
pacf(df$ln_cpue, main="", las=1)
plot_label(lab = "(b)")
mtext(text = "Lag", side = 1, line = 2.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Model autocorrelation assessment:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

par(mfrow=c(2,1), mar=c(3,4.5,1,1)+0.1, oma=c(1,0,0,0))
acf(residuals(mod2), main="", las=1)
plot_label(lab = "(a)", x_prop = 0.92)
pacf(residuals(mod2), main="", las=1)
plot_label(lab = "(b)", x_prop = 0.92)
mtext(text = "Lag", side = 1, line = 2.5)


##############################################################
# Generalized Additive Mixed Models (GAMMs) to account for 
# autocorrelation and cyclic data
##############################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Building a GAMM:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ar_mod = gamm(ln_cpue ~ fInun +  s(water_doy, by=fInun, bs="cc"), data = df,
              correlation = corARMA(form = ~ 1|water_year, p=1, q=1))

# ~ Look at the fitted parameters: phi and theta
summary(ar_mod$lme)$ modelStruct$ corStruct

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ GAMM autocorrelation assessment:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

par(mfrow=c(2,1), mar=c(3,4.5,1,1)+0.1, oma=c(1,0,0,0))
acf(residuals(ar_mod$lme, type = "normalized"), main="", las=1)
plot_label(lab = "(a)", x_prop = 0.92)
pacf(residuals(ar_mod$lme, type = "normalized"), main="", las=1)
plot_label(lab = "(b)", x_prop = 0.92)
mtext(text = "Lag", side = 1, line = 2.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Evaluate GAMM fit:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inun_index_list = list(inun0 = which(df$fInun=="0"),
                       inun1 = which(df$fInun=="1"))

layout(mat = matrix(c(1,1,1,2,2,2,3,3,3,3,4,4), 2, 6, byrow=TRUE))
par(mar=c(4.5,4.5,1.5,1.5)+0.1, oma=rep(0,4))
plot(df$ln_cpue ~ predict(ar_mod$gam), las=1, type="n")
for(i in 1:2){
  points(df$ln_cpue[inun_index_list[[i]]] ~
           predict(ar_mod$lme)[inun_index_list[[i]]], pch=16,
         col=rgb(red = col_mat[1,i],green = col_mat[2,i],
                 blue = col_mat[3,i], alpha=200,
                 maxColorValue = 255))
}
abline(0,1)
plot_label(lab = "(a)")

hist(residuals(ar_mod$lme, type = "normalized"), main="",
     col=rgb(red = col_mat[1,1],green = col_mat[2,1],
             blue = col_mat[3,1], alpha=255,
             maxColorValue = 255), 
     breaks = c(seq(from=-8, to=8, by=0.5)))
hist(residuals(ar_mod$lme, type = "normalized")[inun_index_list[[2]]], add = T,
     col=rgb(red = col_mat[1,2],green = col_mat[2,2],
             blue = col_mat[3,2], alpha=255,
             maxColorValue = 255), 
     breaks = c(seq(from=-8, to=8, by=0.5)))
box(which="plot")
plot_label(lab = "(b)")

plot(residuals(ar_mod$lme, type = "normalized") ~ df$water_doy, las=1, type="n")
for(i in 1:2){
  points(residuals(ar_mod$lme, type = "normalized")[inun_index_list[[i]]] ~
           df$water_doy[inun_index_list[[i]]], pch=16,
         col=rgb(red = col_mat[1,i],green = col_mat[2,i],
                 blue = col_mat[3,i], alpha=200,
                 maxColorValue = 255))
}
abline(h=0)
plot_label(lab = "(c)")

plot(residuals(ar_mod$lme, type = "normalized") ~ df$fInun, las=1, pch=20, 
     xlab="", whisklty=1,
     col=rgb(red = col_mat[1,],green = col_mat[2,],
             blue = col_mat[3,], alpha=255,
             maxColorValue = 255))
abline(h=0)
plot_label(lab = "(d)")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ GAMM Statistical Inference:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Determine the confidence interval to present in your figure and determine the
# t-distribution value of which to multiply the standard error
interval_value = 0.835
t_distribution_probability = 1.0 - (1-interval_value)/2
CI_factor = qt(t_distribution_probability, Inf)


prob_values = c(0, 1)
prob_mat = matrix(data = NA, nrow = 2, ncol = 2, byrow = TRUE,
                  dimnames = list(inundation = c(0,1),
                                  water_doy = paste0("Prob_", prob_values)))
for(i in 1:length(levels(df$fInun))){
  sub = subset(df, df$fInun==levels(df$fInun)[i])
  prob_mat[i,] = round(quantile(x = na.omit(sub$water_doy),
                                probs=prob_values),
                       digits=0)
}

inun_df = data.frame(water_doy = prob_mat[2,1]:prob_mat[2,2],
                     fInun = rep("1", times=
                                   length(prob_mat[2,1]:
                                            prob_mat[2,2])))
dry_df = data.frame(water_doy = prob_mat[1,1]:prob_mat[1,2],
                    fInun = rep("0", times=
                                  length(prob_mat[1,1]:
                                           prob_mat[1,2])))

inun_pred = predict(ar_mod$gam, newdata = inun_df, se.fit = TRUE)
dry_pred = predict(ar_mod$gam, newdata = dry_df, se.fit = TRUE)

inun_fit = inun_pred$fit
inun_upper = inun_fit-CI_factor*inun_pred$se.fit
inun_lower = inun_fit+CI_factor*inun_pred$se.fit

dry_fit = dry_pred$fit
dry_upper = dry_fit-CI_factor*dry_pred$se.fit
dry_lower = dry_fit+CI_factor*dry_pred$se.fit



#~~~~  Stat difference: bottom of the heap:
# https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/

smooth_diff <- function(model, newdata, f1, f2, var, alpha = 0.05,
                        unconditional = FALSE) {
  xp <- predict(model, newdata = newdata, type = 'lpmatrix')
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
  crit <- qt(alpha/2, df.residual(model), lower.tail = TRUE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  data.frame(pair = paste(f1, f2, sep = '-'),
             diff = dif,
             se = se,
             upper = upr,
             lower = lwr)
}
pdat <- expand.grid(water_doy = prob_mat[2,1]:prob_mat[2,2],
                    fInun=c("0", "1"))
comp1 <- smooth_diff(ar_mod$gam, pdat, '0', '1', 'fInun')
comp <- cbind(water_doy = prob_mat[2,1]:prob_mat[2,2],
              comp1)

#Plot it

par(mfrow=c(2,1),mar=c(2.5, 4.5, 2,4.5)+0.1, oma=c(1.5,0,0.5,0))
plot(range(inun_upper, inun_lower, dry_upper,dry_lower) ~
       range(df$water_doy), type = "n", las=1,
     ylab="", xlab="")
plot_label(lab = "(a)")
mtext(text = "CPUE", side = 4, line = 3)
mtext(text = expression(log[e]~"(CPUE)"), side = 2, line = 2.75)
y_labs = c(0.01, 0.02, 0.03, 0.05, 0.08, 0.14, 0.22, 0.37)
axis(side = 4, at = log(y_labs), labels = y_labs, las=1)
polygon(x = c(inun_df$water_doy, max(inun_df$water_doy),
              rev(inun_df$water_doy), inun_df$water_doy[1]),
        y = c(inun_lower, inun_upper[length(inun_upper)], 
              rev(inun_upper), inun_lower[1]), 
        col=rgb(red = col_mat[1,2],green = col_mat[2,2],
                blue = col_mat[3,2], alpha=200,
                maxColorValue = 255),
        border=NA)
lines(inun_fit ~ inun_df$water_doy, lwd=3,
      col="darkblue")
polygon(x = c(dry_df$water_doy, max(dry_df$water_doy),
              rev(dry_df$water_doy), dry_df$water_doy[1]),
        y = c(dry_lower, dry_upper[length(dry_upper)], 
              rev(dry_upper), dry_lower[1]), 
        col=rgb(red = col_mat[1,1],green = col_mat[2,1],
                blue = col_mat[3,1], alpha=180,
                maxColorValue = 255),
        border=NA)
lines(dry_fit ~ dry_df$water_doy, lwd=3,
      col="darkred")
rug(df$water_doy[which(df$fInun=="0")], ticksize = 0.03,
    side = 1, lwd = 0.75,
    col=rgb(red = col_mat[1,1],green = col_mat[2,1],
            blue = col_mat[3,1], alpha=180,
            maxColorValue = 255))
rug(df$water_doy[which(df$fInun=="1")], ticksize = 0.03,
    side = 3, lwd = 0.75,
    col=rgb(red = col_mat[1,2],green = col_mat[2,2],
            blue = col_mat[3,2], alpha=200,
            maxColorValue = 255))
axis(side = 3, at = at_dates, labels = lab_dates)

plot(range(comp$upper, comp$lower) ~ range(df$water_doy),
     type = "n", las=1, ylab="", xlab="")
mtext(text = "Difference in CPUE Trend", side = 2, line = 2.75)
mtext(text = "Day of Water Year", side = 1, line = 2.5)
plot_label(lab = "(b)")
abline(h=0, lty=1, lwd=2)
polygon(x = c(comp$water_doy, max(comp$water_doy), rev(comp$water_doy),
              comp$water_doy[1]),
        y = c(comp$lower, comp$upper[dim(comp)[1]], rev(comp$upper),
              comp$lower[1]), col=gray(0.5,0.5), border=NA)
lines(comp$diff ~ comp$water_doy, lwd=3, col = gray(0.5,1))



