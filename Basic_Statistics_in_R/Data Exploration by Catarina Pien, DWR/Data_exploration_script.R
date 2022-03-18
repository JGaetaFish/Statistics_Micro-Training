### IEP Workshop 2022
### Data Exploration Training


## Covers techniques for exploring a dataset prior to or during analysis
## This is script you can use to follow along during the training.

# Clear environment and load libraries -----------------------
rm(list=ls(all=TRUE))

library(readr) # Read in data
library(dplyr) # Manipulate data
library(tidyr) # Pivot data
library(lubridate) # Datetimes
library(plotly) # Interactive plotting
library(psych) # Correlation plots
library(viridis) # Color palettes
library(ODWGtools) # Tukey and MAD tests
library(leaflet) # Interactive maps
library(ggpubr) # dotplot
source("corvif.R") # From Zuur et al. 2009 for calculating VIF

# Create a theme to make font sizes slightly bigger for plots.----------------
# I have put in an option to rotate the x-axis text when needed (degree = 90 in those cases) (slightly bigger text) 

theme_plots <- function(degree = 0) {
  theme_bw() +
    theme(axis.text.x = element_text(size = 14, angle = degree),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.text = element_text(size = 14),
          strip.text = element_text(size = 11))
}

# Explore data structure -----------------------------------------------

## Read in Yolo Bypass Beach Seine Data -
fishData <- readRDS("ybfmp_fishdata.rds")

## Take a look at the data 
glimpse(fishData)
head(fishData)

# Re-organize data to include and add what you want ----------------------

## Datetime variables, filter years, remove variables
fishSeine <- fishData %>%
  #Calculate CPUE, add datetime variables
  mutate(CPUE = Count/VolumeSeined,
         Date = lubridate::mdy(SampleDate),
         Datetime = paste(Date, SampleTime),
         Datetime = lubridate::ymd_hms(Datetime),
         Year = lubridate::year(Date),
         fYear = factor(Year),
         Month = lubridate::month(Date),
         fMonth = factor(Month),
         WY = ifelse(Month >9, Year + 1, Year),
         fWY = factor(WY))%>%
  # Remove samples in bad condition, filter to a smaller timespan and station list
  filter(GearConditionCode <3,
         Year>2009 & Year<2016)%>%
  # Remove variables not needed
  select(-SampleTime, -SampleDate, -GearID, -GearConditionCode, -MethodCode, -ForkLength, -GeneticallyConfirmed, -GeneticID, -Field_ID_CommonName, -FishSex, -MarkCode, -CWTSample, -StageCode, -Dead, -FishTagID, -Race, -SubstrateCode)


## Modify fish counts to complete counts and calculate daily CPUE/species
fishSeineC <- fishSeine %>%
  # Sum to one count per species per datetime and station (data not formatted this way since some fish have lengths and others do not)
  group_by(StationCode, Datetime, CommonName) %>%
  mutate(Count2 = sum(Count),
         CPUE2 = sum(CPUE),
         Volume2 = sum(VolumeSeined)) %>%
  ungroup() %>%
  select(-Count, -CPUE, -VolumeSeined) %>%
  distinct() %>%
  rename(Count = Count2, CPUE = CPUE2, VolumeSeined = Volume2) %>%
  arrange(CommonName) %>%
  # Add zeroes in
  pivot_wider(names_from = CommonName, values_from = CPUE,values_fill = 0L) %>%
  janitor::clean_names(., case = "lower_camel") %>%
  select(-noCatch) %>%
  pivot_longer(cols = c(americanShad:yellowfinGoby), names_to = "commonName", values_to = "cpue") 

## Filter out stations of interest
fishSeineS <- fishSeineC %>%
  filter(stationCode %in% c("AL2","LI", "BL2", "BL4"))

## Create datasets of interest --------------------
### Define subset native species
nativespecies <- c("chinookSalmon", "splittail", "sacramentoBlackfish", "sacramentoPikeminnow", "sacramentoSucker", "hardhead", "tulePerch", "hitch")

### Filter seine data to species above
fishNative <- filter(fishSeineS, commonName %in% nativespecies)

### Calculate daily native fish CPUE
fishCPUE_n <- fishNative %>%
  group_by(stationCode, date, year, fWy, fMonth) %>%
  summarize(sumCPUE = sum(cpue)) %>%
  ungroup()

### Calculate daily fish CPUE (all fish)
fishCPUE <- fishSeineS %>%
  group_by(stationCode, date, year, fWy, fMonth) %>%
  summarize(sumCPUE = sum(cpue)) %>%
  ungroup()

### Calculate daily native fish CPUE for all stations
fishSeineNatives <- fishSeineC %>%
  filter(commonName %in% nativespecies) %>%
  group_by(wy, date, stationCode, latitude, longitude) %>%
  summarize(sumCPUE = sum(cpue)) %>%
  filter(wy>2011 & wy<2015) %>%
  mutate(fWy = factor(wy),
         fMonth = factor(month(date)))

### Calculate annual native fish CPUE for all stations
fishAnnualNatives <- fishSeineNatives %>%
  group_by(fWy, stationCode, latitude, longitude) %>%
  summarize(sumcpue = sum(sumCPUE),
            n = n())

### Water quality data
fish_WQ <- fishSeineS %>%
  select(date, year, fWy, fMonth, stationCode, waterTemperature, conductivity, do, pH, secchi, turbidity, latitude, longitude) %>%
  distinct() %>%
  mutate(samplingId = 1:nrow(.))

### Water quality in long
WQ_long <- pivot_longer(fish_WQ, cols = waterTemperature:turbidity, names_to = "parameter", values_to = "value") %>%
  mutate(index = 1:nrow(.))

### Join daily fish CPUE with covariates
vars <- left_join(fishCPUE, fish_WQ) %>%
  select(fWy, fMonth, date, stationCode, sumCPUE, waterTemperature, conductivity, turbidity, secchi, do, pH, latitude, longitude)

# Plot data ----------------------------------------------------

## A) Data over time ----------------
# dotchart
dotchart(vars$sumCPUE, bg = "springgreen4")
dotchart(vars$sumCPUE, labels = vars$stationCode, 
         groups = vars$longitude, 
         bg = viridis(4, option = "magma"))
# point plot
(Seine_pointplots <- vars %>%
    ggplot(aes(x= date, y = sumCPUE, color = stationCode)) + 
    geom_point() + 
    geom_vline(xintercept = as.Date("2010-10-01"), linetype = "dashed") + 
    geom_vline(xintercept = as.Date("2011-10-01"), linetype = "dashed") + 
    geom_vline(xintercept = as.Date("2012-10-01"), linetype = "dashed") + 
    geom_vline(xintercept = as.Date("2013-10-01"), linetype = "dashed") + 
    geom_vline(xintercept = as.Date("2014-10-01"), linetype = "dashed") + 
    geom_vline(xintercept = as.Date("2015-10-01"), linetype = "dashed") + 
    viridis::scale_color_viridis(discrete = TRUE, option = "turbo") + 
    scale_x_date(breaks = "6 months") +
    theme_plots(90))


## B) Data over space ----------------
stations <- fishSeineNatives %>%
  select(stationCode, longitude, latitude, sumCPUE) %>%
  group_by(stationCode, longitude, latitude) %>%
  filter(sumCPUE !=0) %>% 
  slice(1) %>% # just picking one value per station for this exercise. slice picks the first value for each station-lon-lat combo.
  mutate(sumCPUE = round(sumCPUE,3)) # just rounding CPUE so it's not a crazy long number when it pops up on leaflet.

# Define your color palette by telling it what colors 
# to use for which variable. You will need one color per level.
pal <- colorFactor(viridis::turbo(13), stations$stationCode)  

### Leaflet interactive map
# Make map
# Make map
stations %>%
  leaflet() %>% # call leaflet.
  addTiles() %>% # this adds the map in the background.
  addCircleMarkers(
    color = ~pal(stationCode),
    stroke = FALSE, # alters whether there is a border to circle
    fillOpacity = 0.9,
    radius = ~sumCPUE*400, # leaflet isn't always great at picking variation in size, here I multiplied by 400 to get something looking appropriate.
    lng = ~longitude, # call your longitude column name
    lat = ~latitude, # call you latitude column name
    label = ~paste(stationCode, " sumCPUE:", sumCPUE,  "Lat:", latitude, "Long:", longitude)) %>% 
  # edit what you want to show up in your label
  addLegend(pal = pal,
            values = ~stationCode,
            position = "bottomright")

### Look at some patterns in the data and note that certain stations have greater CPUE across years while others fluctuate
ggplot(fishAnnualNatives, aes(longitude, latitude, 
                              size = sumcpue, color = stationCode), 
       alpha = 0.5) + 
  geom_jitter() + 
  facet_wrap(~fWy) + 
  scale_colour_viridis(option = "turbo", discrete = TRUE)+ 
  theme_plots(90) +
  theme(legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(ncol = 3, byrow = T),
         size = guide_legend(ncol = 2))

## C) Variability of data and continuous-categorical relationships --------

# We can use boxplots to look at the variability of your continuous data by categorical variables, and also to identify major outliers. 

### By Station -----
(Seine_boxplot <- WQ_long %>%
    ggplot(aes(x= stationCode, y = value, fill = stationCode)) + 
    geom_boxplot() +
    facet_wrap(~parameter, scales = "free") +
    viridis::scale_fill_viridis(discrete = TRUE, option = "turbo") +
    theme_plots(90))

### By Water Year -------
(Seine_boxplots_WY <- WQ_long %>%
    ggplot(aes(x= fWy, y = value, fill = fWy)) + geom_boxplot() +
    facet_wrap(~parameter, scales = "free") +
    viridis::scale_fill_viridis(discrete = TRUE, option = "turbo") +
    theme_plots(90)) 

### By Station and Water Year ------
(Seine_boxplots_month <- WQ_long %>%
    ggplot(aes(x= fWy, y = value, fill = stationCode)) + geom_boxplot() +
    facet_wrap(~parameter, scales = "free") +
    viridis::scale_fill_viridis(discrete = TRUE, option = "turbo") +
    theme_plots(90))

# Fish data
(Fish_boxplots <- fishSeineNatives %>%
    ggplot(aes(x= factor(wy), y = sumCPUE, fill = factor(wy))) + 
    geom_boxplot() +
    viridis::scale_fill_viridis(discrete = TRUE, option = "turbo") +
    theme_plots(90))

## D) Collinearity and relationships between continuous variables------------------

# Scatterplots are useful to look at relationships between continuous variables. These plots are useful:
# 1) to help you find outliers in variables that *should* be correlated
# 2) to help you see relationships that might be interesting to explore
# 3) to help determine whether certain variables may need to be removed for your model 


# Turbidity and Secchi should be correlated
ggplot(vars,aes(turbidity, secchi)) + geom_point() + geom_smooth() + theme_plots()
#DO and Water Temperature should be correlated
ggplot(vars,aes(do, waterTemperature)) + geom_point() + geom_smooth(method = "lm") + theme_plots()
ggplot(vars,aes(do, waterTemperature)) + 
  geom_point() + 
  geom_smooth(method = "lm") + facet_wrap(~stationCode) + theme_plots()

# Pairplots are useful for identifying collinearity (variables that are correlated with each other)
vars_filt <- vars %>% select(sumCPUE, waterTemperature, conductivity, secchi, turbidity, pH, do)
(corr_plot <- pairs.panels(vars_filt, hist.col = "white", cex.cor = 1.5, pch = 21, stars = TRUE))

### Variance Inflation Factor
#If values are below 3, they are acceptable to be included together in a model. 
corvif(vars_filt[,-1])

## E) Distribution of data---------------------
# We can use histograms to look at the distribution of the data, and can also find outliers this way. 

par(mfrow = c(2,2))

ggplot(fishCPUE_n, aes(sumCPUE)) + geom_histogram(binwidth = 0.1) + 
  labs(title = "Native Fish") + theme_plots()

ggplot(fishCPUE, aes(sumCPUE)) + geom_histogram(binwidth = 0.1) + 
  labs(title = "All Fish") + theme_plots()

# freqpoly (multiple categories)
ggplot(fishCPUE, aes(sumCPUE, color = stationCode)) + 
  geom_freqpoly(binwidth = 0.5) + labs(title = "All Fish by Station") + theme_plots()

ggplot(WQ_long, aes(value)) + 
  geom_histogram(color = "navyblue", fill = "steelblue4") + 
  facet_wrap(~parameter, scales = "free") + theme_plots()

# log transform data
par(mfrow = c(1,2))
hist(fish_WQ$turbidity, main = "Turbidity")
hist(log(fish_WQ$turbidity), main = "Log-transformed Turbidity", col = "steelblue4")

## F) How balanced are your data? -------------------

fishCPUE %>%
  count(stationCode)

ggplot(fishAnnualNatives, aes(longitude, latitude, 
                              size = n, fill = stationCode))+ 
  geom_point(alpha = 0.9, shape = 21, color = "black") + 
  facet_wrap(~fWy) + 
  scale_fill_viridis(option = "turbo", discrete = TRUE) + 
  theme_plots(90)  +
  theme(legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(ncol = 3, byrow = T),
         size = guide_legend(ncol = 2))

vars %>%
  count(fWy)

## G) Are there missing data?------------------
# A tile plot can be useful in identifying gaps in data.

samples <- vars %>%
  group_by(stationCode, fWy, fMonth) %>%
  summarize(n = n())
ggplot(samples, aes(x = fWy, y = fMonth, fill = n)) + geom_tile() +
  facet_wrap(~stationCode) + theme_plots(90)

# Identifying Outliers -----------------------------
(Seine_pointplots <- WQ_long %>%
   ggplot(aes(x= date, y = value, 
              color = stationCode, 
              text = paste("index:", index, "samplingID:", samplingId))) + 
   geom_point(size = 0.5) +  
   facet_wrap(~parameter, scales = "free") + 
   viridis::scale_color_viridis(discrete = TRUE, option = "turbo") + 
   theme_plots(90))

# Use plotly to help identify sample IDs with outliers, values, information about point
ggplotly(Seine_pointplots)

# Can identify outliers by filtering
pH_flagID <- WQ_long %>% filter(parameter == "pH" & (value <7.5 | value>9.1))
head(pH_flagID)

# Dotchart
vars_filt$index = row.names(vars_filt)
ggdotchart(data = vars_filt, 
           x= "index", 
           y = "sumCPUE", 
           rotate = T)

## Run outlier tests ------------

outliers <- WQ_long %>%
  group_by(parameter) %>%
  mutate(MAD = ODWGtools::outlier_mad(value),
         Tukey = ODWGtools::outlier_tukey(value))  %>%
  mutate(MAD = ifelse(is.na(MAD), "not outlier", as.character(MAD)),
         Tukey = ifelse(is.na(Tukey), "not outlier", as.character(Tukey)) )

head(outliers %>% select(date, value, MAD, Tukey))

## Visualize outliers ------------------------------
ggplot(outliers) + 
  geom_point(aes(x = date, y = value, color = Tukey, shape = MAD)) + 
  facet_wrap(~parameter, scales = "free") + 
  scale_color_viridis(discrete = TRUE, option = "turbo") + 
  theme_plots(90)

