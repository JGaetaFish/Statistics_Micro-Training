#Ordination introduction
#Tim D. Malinich
#timothy.malinich@wildlife.ca.gov
#Edition 2022

#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#

setwd("C:/Users/tmalinich/Desktop/Ordination_Intro")#Connect to files


#Download the original, if needed, CatchPerStation Flatfile from the STN webpage: https://filelib.wildlife.ca.gov/Public/TownetFallMidwaterTrawl/TNS%20MS%20Access%20Data/TNS%20data/
#Don't forget to examine the Meta-data file when diving into a new dataset!



STN<-read.csv("Supplemental STN_flatfile.csv.",fileEncoding="UTF-8-BOM")#Summer Townet Catch per Station and the Delta Smelt Life Cycle Model Regional Assignments

#always check your dataset
View(STN)#see the dataset in all its glory.
names(STN)#useful when referencing columns
summary(STN$DSLCM.Region)#See the summary of regions

##

library(vegan)  #Statistics
library(ggplot2)#graphing
library(dplyr)  #tools for dataset work
library(tidyr)  #tools for dataset work

#--------------#Datasetprep techniques#-------------------------#
#The next steps I take here are just to make this
#introduction simpler. They are optional for you.

ds<-na.omit(ds)#For this example we'll omit rows with NA's (blanks in the dataset)

ds<-ds[ds$Year>2011,]#lets examine only the last 10 years. Skip if you wish to include past years

#--------------#Ordination techniques#-------------------------#

#Principal Components Analysis (Base R, Stats packages)
#--#prcomp()
#--#princomp()

#Non-Metric Multidimensional Scaling (vegan)
#--#MetaDMS()

#Correspondence Analysis (vegan) 
#--#cca()

#Redundancy Analysis (vegan)
#--#rda()


#------------------------Demonstration--------------------------#
#---------------------------------------------------------------#
#---------------------------------------------------------------#

#Summer Townet Dataset

#many columns for species
#many columns for environmental metrics
names(ds)#Here is the list of possible variables


#--------------------------------------------------------------#
#------question and analysis-----------------------------------#
#--------------------------------------------------------------#

summary(lm(log(ds$Secchi)~ds$DSLCM.Region))#are there differences among my regional variables?

#interesting, but what if we wanted to look at a suite of environmental variables?
# (multivariate)

plot<-ggplot(ds,aes(Temperature.Top,Secchi))+geom_point()
plot#two metrics

plot<-ggplot(ds,aes(Temperature.Top,Secchi,color=DSLCM.Region))+geom_point()
plot#two metrics, and region

#what if I wanted to add Conductivity or more variables?


#Enter ordination techniques (PCA example)

#First lets separate a matrix for our principle Components Analysis

names(ds)#nice to have the names printed out for selecting

ds_m<-ds%>%               # This line sets up a 'pipe' using %>%, meaning that everything following will come from the dataset ds.
  select(Temperature.Top, # This function pulls columns, by name, from the dataset 'ds'
         Secchi,
         Conductivity.Top,
         Depth.Bottom,
         Microcystis,
         Turbidity.Top,
         Weather,
         Waves) 

#This matrix, ds_m, we have created has 8 columns, which means our PCA should produce 8 PCs.


#-------------------------------------------------------------#
#-------------------------------------------------------------#
#---------------------The PCA---------------------------------#

pca_ds<-prcomp(ds_m,center=TRUE,scale.=TRUE)#run the PCA

names(pca_ds)#check out the components of pca_ds

summary(pca_ds)#examine the summary of Principal components
#Look at the cumulative proportion of variance. 1 PC for each column (col=8).
#whenever graphing these PC's, display the % variance explained.

PCscores<-as.data.frame(pca_ds$x)         #scores, save them here for plotting
rotation<-as.data.frame(pca_ds$rotation)  #is the relationship between the original values and PC scores
vectors<-rownames(rotation)               #save here for printing in plot later


#----------------------------------------------------------------------------------------#
#----------plotting principal components--------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#Plot components, examine eigenvectors to see what patterns emerge. 



#Here is an example of the first two PCs and eigenvectors

figure1<-ggplot(data=PCscores,aes(PC1,PC2))+
  geom_segment(data=rotation,aes(x = 0,xend=PC1*10, y =0,yend= PC2*10),
               arrow = arrow(),size=1)+
  geom_text(data=rotation,aes(PC1*10+1,PC2*10+1),label=vectors,size=4)+
  
  xlab("PC 1 (37%)")+ylab("PC 2 (15%)")+
  
  theme_bw()+
  theme(axis.ticks=element_line(color="black"),
        axis.title.x=element_text(size=16, color="Black"),
        axis.title.y=element_text(size=16, color="Black"),
        axis.text.x=element_text(size=14, color = "Black",angle=90),
        axis.text.y=element_text(size=14, color="Black"))
figure1



############different axes

figure2<-ggplot(data=PCscores,aes(PC2,PC3))+
  geom_segment(data=rotation,aes(x = 0,xend=PC2*10, y =0,yend= PC3*10),
               arrow = arrow(),size=1)+
  geom_text(data=rotation,aes(PC2*10+1,PC3*10+1),label=vectors,size=4)+
  
  xlab("PC 2 (15%)")+ylab("PC 3 (11%)")+
  theme_bw()+
  theme(axis.ticks=element_line(color="black"),
        axis.title.x=element_text(size=16, color="Black"),
        axis.title.y=element_text(size=16, color="Black"),
        axis.text.x=element_text(size=14, color = "Black",angle=90),
        axis.text.y=element_text(size=14, color="Black"))
figure2


################################################
########lets add our points
figure3<-ggplot()+
  geom_point(data=PCscores,aes(PC1,PC2),alpha=.3)+
  geom_segment(data=rotation,aes(x = 0,xend=PC1*10, y =0,yend= PC2*10),
               arrow = arrow(),size=1.0)+
  geom_text(data=rotation,aes(PC1*10+1,PC2*10+1),label=vectors,size=4)+
  xlab("PC 1 (37%)")+ylab("PC 2 (15%)")+
  theme_bw()+
  theme(axis.ticks=element_line(color="black"),
        axis.title.x=element_text(size=16, color="Black"),
        axis.title.y=element_text(size=16, color="Black"),
        axis.text.x=element_text(size=14, color = "Black",angle=90),
        axis.text.y=element_text(size=14, color="Black"))
figure3


###################################################

figure4<-ggplot()+
  geom_point(data=PCscores,aes(PC1,PC2,color=ds$DSLCM.Region),alpha=.2)+
  geom_segment(data=rotation,aes(x = 0,xend=PC1*10, y =0,yend= PC2*10),
               arrow = arrow(),size=1.0)+
  geom_text(data=rotation,aes(PC1*10+1,PC2*10+1),label=vectors,size=4)+
  xlab("PC 1 (37%)")+ylab("PC 2 (15%)")+
  theme_bw()+
  theme(axis.ticks=element_line(color="black"),
        axis.title.x=element_text(size=16, color="Black"),
        axis.title.y=element_text(size=16, color="Black"),
        axis.text.x=element_text(size=14, color = "Black",angle=90),
        axis.text.y=element_text(size=14, color="Black"))
figure4


#################use different Principal components

figure5<-ggplot()+
  geom_point(data=PCscores,aes(PC1,PC5,color=ds$DSLCM.Region),alpha=.2)+
  geom_segment(data=rotation,aes(x = 0,xend=PC1*10, y =0,yend= PC5*10),
               arrow = arrow(),size=1.0)+
  geom_text(data=rotation,aes(PC1*10+1,PC5*10+1),label=vectors,size=4)+
  xlab("PC 1 (37%)")+ylab("PC 5 (9%)")+
  theme_bw()+
  theme(axis.ticks=element_line(color="black"),
        axis.title.x=element_text(size=16, color="Black"),
        axis.title.y=element_text(size=16, color="Black"),
        axis.text.x=element_text(size=14, color = "Black",angle=90),
        axis.text.y=element_text(size=14, color="Black"))
figure5


##############continue examining different components

figure6<-ggplot()+
  geom_point(data=PCscores,aes(PC1,PC8,color=ds$DSLCM.Region))+
  geom_segment(data=rotation,aes(x = 0,xend=PC1*10, y =0,yend= PC8*10),
               arrow = arrow(),size=1.5)+
  geom_text(data=rotation,aes(PC1*10+1,PC8*10+1),label=vectors,size=4)+
  xlab("PC 1 (37%)")+ylab("PC 8 (5%)")+
  theme_bw()+
  theme(axis.ticks=element_line(color="black"),
        axis.title.x=element_text(size=16, color="Black"),
        axis.title.y=element_text(size=16, color="Black"),
        axis.text.x=element_text(size=14, color = "Black",angle=90),
        axis.text.y=element_text(size=14, color="Black"))
figure6


##############use ellipses to distinguish regions

figure7<-ggplot()+
  geom_point(data=PCscores,aes(PC1,PC5,color=ds$DSLCM.Region),size=1)+
  stat_ellipse(data=PCscores,aes(PC1,PC5,color=ds$DSLCM.Region),size=1.5)+#ellipse plotted based on multivariate t distribution 
  geom_segment(data=rotation,aes(x = 0,xend=PC1*10, y =0,yend= PC5*10),
               arrow = arrow(),size=1.0)+
  geom_text(data=rotation,aes(PC1*10+1,PC5*10+1),label=vectors,size=4)+
  scale_color_discrete(name="")+
  xlab("PC 1 (37%)")+ylab("PC 5 (9%)")+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(axis.ticks=element_line(color="black"),
        axis.title.x=element_text(size=16, color="Black"),
        axis.title.y=element_text(size=16, color="Black"),
        axis.text.x=element_text(size=14, color = "Black",angle=90),
        axis.text.y=element_text(size=14, color="Black"))
figure7


######################################
############Plenty of other R packages that can help with plotting your PCA results
library(ggfortify)

autoplot(pca_ds,loadings=TRUE,loadings.label=TRUE)

plot(pca_ds$x)#base R package

######################################Do I have good groups?

summary(lm(PCscores$PC1~ds$DSLCM.Region))
summary(manova(as.matrix(PCscores[,1:3])~ds$DSLCM.Region))

#to identify better regions, I suggest RDA or DFA. 

######################################Can use principle components in further analyses

summary(lm(Age.0.Striped.Bass~Conductivity.Top+Turbidity.Top+Waves+Weather+Microcystis+Temperature.Top+Depth.Bottom+Secchi+Tridentiger.spp,data=ds))
#R2 low, but most of this value may be due to the increased number of variables


summary(lm(ds$Age.0.Striped.Bass~PCscores$PC1+ds$Tridentiger.spp))#Account for environmental metrics, but examine the relationship between invasive gobies and SB
#note the lower R2 value

summary(lm(ds$Age.0.Striped.Bass~PCscores$PC1+PCscores$PC2+ds$Tridentiger.spp))#Account for environmental metrics, but examine the relationship between invasive gobies and longfin
#note the still low but improved R2 value

#You have reached the end of this R Script
#You are rewarded with this fish
#
#........................................
#....................<><.................
#........................................