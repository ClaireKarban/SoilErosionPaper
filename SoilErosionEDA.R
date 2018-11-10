#Using dplyr for data manipulation and ggplot for plotting data
#In this script, I will read in data and explore distributions and relationships

#Libraries:
library(ggplot2)
library(dplyr)
library(lme4)   
library(arm)
library(gridExtra)

erosion_dat <- read.csv("SoilErosion.csv", as.is=TRUE) #read in data

head(erosion_dat) #inspect the data

BSNE_dat <- erosion_dat %>%
    filter(BSNE_gM2day != "#N/A") #create new wind erosion data set with the the NA values removed

str(BSNE_dat)


##### Edit variables to make sure they are in the correct format
BSNE_dat$BSNE_gM2day <- as.numeric(BSNE_dat$BSNE_gM2day)
BSNE_dat$Chla <- as.numeric(BSNE_dat$Chla)
BSNE_dat$Site <- as.factor(BSNE_dat$Site)
BSNE_dat$Treatment <- as.factor(BSNE_dat$Treatment)
BSNE_dat$Seeding <- as.factor(BSNE_dat$Seeding)
BSNE_dat$Treatment_Year <- as.factor(BSNE_dat$Treatment_Year)

str(SMBSNE_dat) #Check structure to make sure the variable all look good


##Graphical views of the data - Wind Erosion
ggplot(data=BSNE_dat) +
  geom_histogram(mapping = aes(x=BSNE_gM2day),bins=50)
ggplot(data=BSNE_dat) +
  geom_histogram(mapping = aes(x=BSNE_gM2day, col=Site), bins=50, alpha=0.5) #This is kind of a shitty graph, but it's clear that the values are different.
#I think this is justification for running two separate models.
ggplot(data=BSNE_dat) +
  geom_histogram(mapping = aes(x=BSNE_gM2day), bins=50) +
  facet_wrap(facets = ~Site)

#Data is clearly non-normally distributed. A log transformation should be appropriate 

ggplot(data=BSNE_dat) +
  geom_histogram(mapping = aes(x=BSNE_gM2day), bins=50) +
  facet_wrap(facets = ~Site)

#Add column of log-trasformed BSNE, adding 0.1 to the zeros before transformation
BSNE_dat <- mutate(BSNE_dat,log_BSNE_gM2day=log(ifelse(BSNE_gM2day==0,0.1,BSNE_gM2day))) 

SMBSNE_dat <- subset(BSNE_dat, Site=="SM") #Create new data sets by site so I run sepearate models by site
WMBSNE_dat <- subset(BSNE_dat, Site=="WM")


#plot log-transformed data 
SMBSNE_dat %>%
  group_by(Site) %>%
  ggplot() +
  geom_histogram(mapping = aes(x=log_BSNE_gM2day), bins=50)   #Looks more normal. There's still an issue with lots of zeros

ggplot(data=WMBSNE_dat) +
  geom_histogram(mapping = aes(x=log_BSNE_gM2day), bins=50) #Some values are negative? That's weird


ggplot(data=SMBSNE_dat) +
  geom_histogram(mapping = aes(x=log_BSNE_gM2day, fill=Treatment_Year), position="identity",
                 bins=50, alpha=0.5) #This result is SO weird to me. There was so little erosion in year one for some samples. 

#Group by treatment to see if this is a result of seeding disturbance (i.e. are those high yr 1 disturbances because of seeding)
ggplot(data=SMBSNE_dat) +
  geom_histogram(mapping = aes(x=log_BSNE_gM2day, fill=Treatment_Year), position="identity",
                 bins=50, alpha=0.5) +
  facet_wrap(facets=~Seeding) #Nope, the seeded treatments seem to have the lower erosion. 

ggplot(data=SMBSNE_dat) +
  geom_histogram(mapping = aes(x=log_BSNE_gM2day, fill=Treatment_Year), position="identity",
                 bins=50, alpha=0.5) +
  facet_wrap(facets=~Treatment_Year) 
#Wray
ggplot(data=WMBSNE_dat) +
  geom_histogram(mapping = aes(x=log_BSNE_gM2day, fill=Treatment_Year), position="identity",
                 bins=50, alpha=0.5) #This looks more like what I'd expect to see

ggplot(data=WMBSNE_dat) +
  geom_histogram(mapping = aes(x=log_BSNE_gM2day, fill=Treatment_Year), position="identity",
                 bins=50, alpha=0.5) +
  facet_wrap(facets = ~Seeding)

ggplot(data=WMBSNE_dat) +
  geom_histogram(mapping = aes(x=log_BSNE_gM2day, fill=Treatment_Year), position="identity",
                 bins=50, alpha=0.5) +
  facet_wrap(facets = ~Treatment_Year)

###Silt Erosion
erosion_dat<- mutate(erosion_dat,log_Silt_Erosion=log(ifelse(Silt_Erosion==0,0.1,Silt_Erosion))) 
SMSilt_dat <- subset(erosion_dat, Site=="SM") #Create new data sets by site so I run sepearate models by site
WMSilt_dat <- subset(erosion_dat, Site=="WM")

#plot log-transformed data 
SMSilt_dat %>%
  group_by(Site) %>%
  ggplot() +
  geom_histogram(mapping = aes(x=Silt_Erosion), bins=50)

SMSilt_dat %>%
  group_by(Site) %>%
  ggplot() +
  geom_histogram(mapping = aes(x=log_Silt_Erosion), bins=50)   #Looks more normal but NOT normal

ggplot(data=WMSilt_dat) +
  geom_histogram(mapping = aes(x=Silt_Erosion), bins=50) 

ggplot(data=WMSilt_dat) +
  geom_histogram(mapping = aes(x=log_Silt_Erosion), bins=50) #Not looking great

#Now look at split up by year
ggplot(data=WMSilt_dat) +
  geom_histogram(mapping = aes(x=Silt_Erosion, fill=Treatment_Year), position="identity",
                 bins=50, alpha=0.5) +
  facet_wrap(facets = ~Treatment_Year)

ggplot(data=SMSilt_dat) +
  geom_histogram(mapping = aes(x=Silt_Erosion, fill=Treatment_Year), position="identity",
                 bins=50, alpha=0.5) +
  facet_wrap(facets = ~Treatment_Year) #There does appear to be slightly more erosion in year 

#Before I split the data up by site, are there any differences between Year versus Treatment_year?
ggplot(data=erosion_dat) +
  geom_point(mapping = aes(x=Silt_Erosion, fill=Treatment_Year), position="identity",
                 bins=50, alpha=0.5) +
  facet_wrap(facets = ~Year) 


ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Year, y=BSNE_gM2day, col=Site), position="jitter") +
  labs(x="Year",y="Wind Erosion (g/M2/day)")
  #This is a nice graph showing that Treatment year and site look much more influential on wind erosion than Year.

ggplot(data=erosion_dat) +
  geom_point(mapping = aes(x=Year, y=Silt_Erosion, col=Site), position="jitter") +
  labs(x="Year",y="Water Erosion (g/M2)")
#Notes: Shay has more wind and water erosion than Wray, so the effect of year 1 versus year 2 is more dramatic at Shay.

###Finally, graph predictor variables with response variable
#Wind Erosion
ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Median_SoilStability, y=log(BSNE_gM2day+.1), col=Site), position="jitter") #Not a tight relationship, but looks negative and linear

ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=VegCover, y=log(BSNE_gM2day+.1), col=Site), position="jitter") #This is a better relationship; negative and linear

ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Chla, y=log(BSNE_gM2day+.1), col=Site), position="jitter") #There looks like NO relationship here for WM, and maybe a weakly negative one for SM
ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Chla, y=BSNE_gM2day, col=Site), position="jitter") #This non-log doesn't help much.


ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=P_Grasses._Sedges, y=log(BSNE_gM2day+.1), col=Site), position="jitter") #Not a strong relationship; maybe stronger for WM than SM
ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=P_Grasses._Sedges, y=BSNE_gM2day, col=Site), position="jitter") #Without the log, it looks like more of a relationship, especially for SM!


ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=A_Grasses, y=log(BSNE_gM2day+.1), col=Site), position="jitter") #Might be interesting for SM, but WM just doesn't have any grasses growing

ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=A_Grasses, y=BSNE_gM2day, col=Site), position="jitter")


ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Chla, y=Treatment_Year, col=Site), position="jitter") #Chla should be correlation to years of recovery, but it isn't closely..
ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Chla, y=LightCrust, col=Site), position="jitter") #Okay, this relationship is more as I'd expect; positive linear relationship between light crust and chl a 
ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Chla, y=DarkCrust, col=Site), position="jitter") #Less of a correlation with dark crust
ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Lichen, y=DarkCrust, col=Site), position="jitter") #Dark crust also isn't as correlated to moss or lichen as I'd expect

#Exploring some more ground cover variables
ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Litter, y=log(BSNE_gM2day+.1), col=Site), position="jitter") #This is surprising! This relationship looks pretty good/negative for WM but not for SM at all!
ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Litter, y=log(BSNE_gM2day+.1), col=Treatment), position="jitter") #Lots of litter for B, followed by P and C. L had no litter and high erosion!

ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=LightCrust, y=log(BSNE_gM2day+.1), col=Site), position="jitter")#This looks important, especially for SM!

ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=DarkCrust, y=log(BSNE_gM2day+.1), col=Site), position="jitter") #Less important but still looks like there's a negative relationship

ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Soil, y=log(BSNE_gM2day+.1), col=Site), position="jitter") #Positive relationship, especially for SM!
ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Soil, y=log(BSNE_gM2day+.1), col=Treatment), position="jitter") #Much more bare soil for L and P than B and C! Looks positive!


ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Wood, y=log(BSNE_gM2day+.1), col=Site), position="jitter") #Way more wood at WM than SM. possible negative relationship for WM

ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Burnt_RootShrub, y=log(BSNE_gM2day+.1), col=Site), position="jitter")#Not much here; lot's of zeros...
ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Burnt_RootShrub, y=BSNE_gM2day, col=Site), position="jitter")#Not much here; lot's of zeros...


ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Burnt_PileSoil, y=log(BSNE_gM2day+.1), col=Site), position="jitter") 
ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Burnt_PileSoil, y=log(BSNE_gM2day+.1), col=Treatment), position="jitter") #This is kind of hard to interpret, and doesn't look super significant, but it does break up by Treatment very well.
ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Burnt_PileSoil, y=BSNE_gM2day, col=Site), position="jitter") 
ggplot(data=BSNE_dat) +
  geom_point(mapping = aes(x=Burnt_PileSoil, y=BSNE_gM2day, col=Treatment), position="jitter") 
