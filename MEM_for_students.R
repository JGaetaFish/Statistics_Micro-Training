##############################################################
# Linear Mixed Effects Modeling
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

library(lme4)
library(effects)

# Custom function to add text to the same position on every panel:
plot_label = function(lab="(a)", x_prop=0.08, y_prop=0.92,
                      font_type=2, fcex=1.15, usr=par('usr')){
  x_val = usr[1]+(usr[2]-usr[1])*x_prop
  y_val = usr[3]+(usr[4]-usr[3])*y_prop
  text(x = x_val, y = y_val, labels = lab, font = font_type, cex=fcex)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Load and explore the data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat3 = read.csv(file = "USGS_SFBWQ_3stations.csv", header=T)

head(dat3)

# set up to for loop through the stations
plot_sta = unique(dat3$Station_Number)
col_mat = matrix(data = c(140,81,10, 1,102,94, 118,42,131),
                 nrow = 3, ncol = 3, byrow = TRUE)

par(mfrow=c(1,2), oma=c(0,1.1,0,0), mar=c(3.5,2,0.75,0.75)+0.1)
plot(dat3$Discrete_Oxygen ~ dat3$Temperature, type="p", las=1, 
     ylab="", xlab = "",
     pch=20, col=gray(0.25,0.5))
plot_label(lab = "(a)", x_prop = 0.07)
mtext(text = "Dissolved Oxygen", side = 2, line = 2)
mtext(text = "Water Temperature", side = 1, line = 2.25)

plot(dat3$Discrete_Oxygen ~ dat3$Temperature, type="n", las=1, 
     ylab="", xlab = "")
mtext(text = "Water Temperature", side = 1, line = 2.25)
for(i in 1:length(plot_sta)){
  sub = subset(dat3, dat3$Station_Number == plot_sta[i])
  points(sub$Discrete_Oxygen ~ sub$Temperature, pch=20,
         col=rgb(col_mat[i,1],col_mat[i,2],col_mat[i,3],
                 maxColorValue = 255, alpha = 200))
}
plot_label(lab = "(b)", x_prop = 0.07)
box(which="outer")

##############################################################
# Random Intercept Mixed Effects Regression
##############################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Build the linear model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mod = lm(Discrete_Oxygen ~ Temperature, data = dat3)
summary(mod)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Build the random intercept model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# your grouping variable must be a factor:
dat3$fStation = factor(dat3$Station_Number)

mem = lmer(Discrete_Oxygen ~ Temperature + (1|fStation),
           data = dat3, REML = TRUE)
summary(mem)

# extract model variance terms:
sigma_mem = data.frame(VarCorr(mem))[,5]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Plot the data and both models:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create a vector that spans the observed range of temperatures for simulations
temp_seq = seq(min(dat3$Temperature), max(dat3$Temperature),
               length.out=100)

# Use the effects package to calculate a 95% CI around the grand mean model

mem_eff=Effect(focal.predictors = "Temperature", mod=mem, 
               xlevels=list(Temperature=temp_seq))
lm_eff=Effect(focal.predictors = "Temperature", mod=mod, 
              xlevels=list(Temperature=temp_seq))

# Use the random effect coefficients to generate random effect predictions 
# across the observed range of temperatures
y_15 = coef(mem)$fStation["15",1] + coef(mem)$fStation["15",2]*temp_seq
y_36 = coef(mem)$fStation["36",1] + coef(mem)$fStation["36",2]*temp_seq
y_649 = coef(mem)$fStation["649",1] + coef(mem)$fStation["649",2]*temp_seq

# plot it
par(mfrow=c(1,1), mar=c(4, 4.5, 1, 1)+0.1, oma=rep(0,4))
plot(dat3$Discrete_Oxygen ~ dat3$Temperature, type="n", las=1, 
     ylab="Dissolved Oxygen", xlab = "Water Temperature")
for(i in 1:length(plot_sta)){
  sub = subset(dat3, dat3$Station_Number == plot_sta[i])
  points(sub$Discrete_Oxygen ~ sub$Temperature, pch=20,
         col=rgb(col_mat[i,1],col_mat[i,2],col_mat[i,3],
                 maxColorValue = 255, alpha = 200))
}
lines(temp_seq, y_36, type="l", lwd=1.5,
      col=rgb(col_mat[1,1],col_mat[1,2],col_mat[1,3],
              maxColorValue = 255, alpha = 225))
lines(temp_seq, y_15, type="l", lwd=1.5,
      col=rgb(col_mat[2,1],col_mat[2,2],col_mat[2,3],
              maxColorValue = 255, alpha = 225))
lines(temp_seq, y_649, type="l", lwd=1.5,
      col=rgb(col_mat[3,1],col_mat[3,2],col_mat[3,3],
              maxColorValue = 255, alpha = 225))

polygon(x = c(temp_seq, temp_seq[100], rev(temp_seq), temp_seq[1]),
        y=c(mem_eff$lower, mem_eff$upper[100], 
            rev(mem_eff$upper), mem_eff$lower[1]),
        col=gray(0.25,0.25), 
        border=NA)

polygon(x = c(temp_seq, temp_seq[100], rev(temp_seq), temp_seq[1]),
        y=c(lm_eff$lower, lm_eff$upper[100], 
            rev(lm_eff$upper), lm_eff$lower[1]),
        col=rgb(0.7,0.1,0.1,maxColorValue = 1,alpha = 0.5), 
        border=NA)
lines(mem_eff$fit~temp_seq, lwd=3.5, col="black")
lines(lm_eff$fit~temp_seq, lwd=3.5, lty=2,
      col=rgb(0.5,0.1,0.1,maxColorValue = 1,alpha = 1))

box(which="outer")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Evaluate and compare model fit:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#generate mixed effects model predictions and residuals
new_dat = data.frame(Temperature = dat3$Temperature,
                     fStation = dat3$fStation)
mem_pred=predict(object = mem, newdata = new_dat)
mem_resid = c(dat3$Discrete_Oxygen - mem_pred)


lm_hist = hist(residuals(mod), breaks=10, plot = FALSE)
lmer_hist = hist(mem_resid, breaks=10, plot = FALSE)
freq_max = max(c(lm_hist$counts, lmer_hist$counts))*1.05

#quartz()
par(mfrow=c(3,2), mar=c(4, 2.75, 1, 1)+0.1, oma=c(0,2,2.5,0))
plot(dat3$Discrete_Oxygen ~ fitted(mod), type="n", las=1, 
     cex.axis=1.1, ylab="", xlab = "")
plot_label(lab = "(a)")
mtext(text = "Observed", side = 2, line = 2.75)
mtext(text = "Predicted", side = 1, line = 2.5)
for(i in 1:length(plot_sta)){
  index = which(dat3$Station_Number == plot_sta[i])
  points(dat3$Discrete_Oxygen[index] ~ fitted(mod)[index], pch=20,
         col=rgb(col_mat[i,1],col_mat[i,2],col_mat[i,3],
                 maxColorValue = 255, alpha = 200))
}
abline(0,1)
mtext(text = "Simple Linear Regression", side = 3, line = 1, font=2)

plot(dat3$Discrete_Oxygen ~ mem_pred, type="n", las=1, 
     cex.axis=1.1, ylab="", xlab = "")
plot_label(lab = "(b)")
mtext(text = "Predicted", side = 1, line = 2.5)
for(i in 1:length(plot_sta)){
  index = which(dat3$Station_Number == plot_sta[i])
  points(dat3$Discrete_Oxygen[index] ~ mem_pred[index], pch=20,
         col=rgb(col_mat[i,1],col_mat[i,2],col_mat[i,3],
                 maxColorValue = 255, alpha = 200))
}
abline(0,1)
mtext(text = "Random Intercept\nMixed Effects Regression", side = 3, line = 0.5, font=2)

hist(residuals(mod), breaks=10, main="", xlab="", ylab="",cex.axis=1.1,
     xlim=c(-3.35, 3.35), las=1, yaxs="i", ylim=c(0,freq_max))
box(which="plot")
plot_label(lab = "(c)")
mtext(text = "Frequency", side = 2, line = 2.75)
mtext(text = "Residual", side = 1, line = 2.5)

hist(mem_resid, breaks=10, main="", 
     xlab="", ylab="",cex.axis=1.1,
     xlim=c(-3.35, 3.35), las=1, yaxs="i", ylim=c(0,freq_max))
box(which="plot")
plot_label(lab = "(d)")
mtext(text = "Residual", side = 1, line = 2.5)

#~~~~~~~~~~~~~~~~~~~~~~
plot(dat3$Temperature , residuals(mod), type="n", las=1, 
     ylab="", xlab = "",cex.axis=1.1,
     ylim=range(c(residuals(mod), mem_resid)))
plot_label(lab = "(e)")
mtext(text = "Residual", side = 2, line = 2.75)
mtext(text = "Water Temperature", side = 1, line = 2.5)
for(i in 1:length(plot_sta)){
  index = which(dat3$Station_Number == plot_sta[i])
  points(dat3$Temperature[index] , residuals(mod)[index], pch=20,
         col=rgb(col_mat[i,1],col_mat[i,2],col_mat[i,3],
                 maxColorValue = 255, alpha = 200))
}
abline(h=0)

plot(dat3$Temperature , mem_resid, 
     type="n", las=1, cex.axis=1.1,
     ylim=range(c(residuals(mod), mem_resid)),
     ylab="", xlab = "")
plot_label(lab = "(f)")
mtext(text = "Water Temperature", side = 1, line = 2.5)
for(i in 1:length(plot_sta)){
  index = which(dat3$Station_Number == plot_sta[i])
  points(dat3$Temperature[index] , 
         mem_resid[index], pch=20,
         col=rgb(col_mat[i,1],col_mat[i,2],col_mat[i,3],
                 maxColorValue = 255, alpha = 200))
}
abline(h=0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Visualize model uncertainty:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

varCor = as.data.frame(VarCorr(mem))
varCor = varCor[is.na(varCor$var2),]

int_rnorm = rnorm(n = 100000, mean = fixef(mem)[1],
                  sd = varCor$sdcor[which(varCor$grp=="fStation" & 
                                            varCor$var1=="(Intercept)")])

int_rnorm_den = density(int_rnorm, cut = 0, adjust=1)
int_den = density(coef(mem)$fStation[,1], cut=0, adjust=0.25)

int_rnorm_den_df = data.frame(x = int_rnorm_den$x,
                              y_nu = int_rnorm_den$y/max(int_rnorm_den$y))
int_den_df = data.frame(x = int_den$x,
                        y_nu = int_den$y/max(int_den$y))

resid_rnorm = rnorm(n = 100000, mean = 0,
                    sd = varCor$sdcor[which(varCor$grp=="Residual")])

resid_mem = residuals(mem)
resid_mem_15 = resid_mem[which(dat3$fStation=="15")]
resid_mem_36 = resid_mem[which(dat3$fStation=="36")]
resid_mem_649 = resid_mem[which(dat3$fStation=="649")]

resid_rnorm_den = density(resid_rnorm, cut = 0, adjust=1)
resid_15_den = density(resid_mem_15, cut=0, adjust=0.75)
resid_36_den = density(resid_mem_36, cut=0, adjust=0.75)
resid_649_den = density(resid_mem_649, cut=0, adjust=0.75)

resid_rnorm_den_df = data.frame(x = resid_rnorm_den$x,
                                y_nu = resid_rnorm_den$y/max(resid_rnorm_den$y))
resid_15_den_df = data.frame(x = resid_15_den$x,
                             y_nu = resid_15_den$y/max(resid_15_den$y))
resid_36_den_df = data.frame(x = resid_36_den$x,
                             y_nu = resid_36_den$y/max(resid_36_den$y))
resid_649_den_df = data.frame(x = resid_649_den$x,
                              y_nu = resid_649_den$y/max(resid_649_den$y))


par(mfrow=c(1,2), mar=c(3.5,3,2,1.5)+0.1, oma=c(0,1,1,0))
plot(int_rnorm_den_df, ylim=c(0,1.05), las=1, type="n",
     yaxs="i", xlab="", ylab="", cex.axis=1)
mtext(text = "Relative Frequency", side = 2, line = 3, cex=1)
mtext(text = "Random Intercept", side = 1, line = 2, cex=1)
polygon(x = c(int_rnorm_den_df$x, int_rnorm_den_df$x[512],
              int_rnorm_den_df$x[1], int_rnorm_den_df$x[1]),
        y = c(int_rnorm_den_df$y_nu, 0, 0,
              int_rnorm_den_df$y_nu[1]), border=NA, col=gray(0.25,0.4))
polygon(x = c(int_den_df$x, int_den_df$x[512],
              int_den_df$x[1], int_den_df$x[1]),
        y = c(int_den_df$y_nu, 0, 0,
              int_den_df$y_nu[1]), border=NA,
        col=rgb(0.68,0.17,0.17,maxColorValue = 1,alpha = 0.66))
plot_label(lab = "(a)")

plot(resid_rnorm_den_df, ylim=c(0,1.05), las=1, type="n",
     yaxs="i", xlab="", ylab="", cex.axis=1)
mtext(text = "Residual", side = 1, line = 2, cex=1)
polygon(x = c(resid_rnorm_den_df$x, resid_rnorm_den_df$x[512],
              resid_rnorm_den_df$x[1], resid_rnorm_den_df$x[1]),
        y = c(resid_rnorm_den_df$y_nu, 0, 0,
              resid_rnorm_den_df$y_nu[1]), border=NA, col=gray(0.25,0.4))
polygon(x = c(resid_36_den_df$x, resid_36_den_df$x[512],
              resid_36_den_df$x[1], resid_36_den_df$x[1]),
        y = c(resid_36_den_df$y_nu, 0, 0,
              resid_36_den_df$y_nu[1]), border=NA,
        col=rgb(col_mat[1,1],col_mat[1,2],col_mat[1,3],
                maxColorValue = 255, alpha = 180))

polygon(x = c(resid_15_den_df$x, resid_15_den_df$x[512],
              resid_15_den_df$x[1], resid_15_den_df$x[1]),
        y = c(resid_15_den_df$y_nu, 0, 0,
              resid_15_den_df$y_nu[1]), border=NA,
        col=rgb(col_mat[2,1],col_mat[2,2],col_mat[2,3],
                maxColorValue = 255, alpha = 180))

polygon(x = c(resid_649_den_df$x, resid_649_den_df$x[512],
              resid_649_den_df$x[1], resid_649_den_df$x[1]),
        y = c(resid_649_den_df$y_nu, 0, 0,
              resid_649_den_df$y_nu[1]), border=NA,
        col=rgb(col_mat[3,1],col_mat[3,2],col_mat[3,3],
                maxColorValue = 255, alpha = 180))
plot_label(lab = "(b)")

par(fig = c(0, 1, 0, 1), oma = c(0, 4, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", xpd=NA, legend = c("Underlying Distribution",  "Random Intercept"), pch=22,
       pt.bg=c(gray(0.3,0.3),
               rgb(0.68,0.17,0.17,maxColorValue = 1,alpha = 0.66)),
       horiz=T, pt.cex = 2, bty="n", inset=-0.01, cex=1)
legend("top", xpd=NA, pch=22, cex=1,
       legend = c("Calaveras Point", "Point San Pablo","Sacramento River"),
       pt.bg=c(rgb(col_mat[1,1],col_mat[1,2],col_mat[1,3],
                   maxColorValue = 255, alpha = 225),
               rgb(col_mat[2,1],col_mat[2,2],col_mat[2,3],
                   maxColorValue = 255, alpha = 225),
               rgb(col_mat[3,1],col_mat[3,2],col_mat[3,3],
                   maxColorValue = 255, alpha = 225)),
       horiz=T, pt.cex = 3, inset=c(0.06), bty="n")
box(which="outer")


##############################################################
# Random Slope Mixed Effects Regression
##############################################################

mem2 = lmer(Discrete_Oxygen ~ Temperature + (Temperature-1|fStation),
           data = dat3, REML = TRUE)
summary(mem2)


##############################################################
# Random Intercept and Slope Mixed Effects Regression
##############################################################

mem3 = lmer(Discrete_Oxygen ~ Temperature + (1+Temperature|fStation),
           data = dat3, REML = TRUE)
summary(mem3)
