##In this script, I will fit partial-pooling multilevel models with lmer() and Bayesian approaches
##I am exploring predictors of wind and water erosion at two sites on the Colorado Plateau
##Sites will be analyzed separately and together.

#Libraries:
library(ggplot2)
library(dplyr)
library(lme4)    



##Tidy up data
erosion_dat <- read.csv("SoilErosion.csv", as.is=TRUE) #read in data

head(erosion_dat) #inspect the data
str(erosion_dat)

erosion_dat$Site <- as.factor(erosion_dat$Site)
erosion_dat$Treatment <- as.factor(erosion_dat$Treatment)
erosion_dat$Chla <- as.numeric(erosion_dat$Chla)
erosion_dat$Treatment <- factor(erosion_dat$Treatment, levels = c("C","B","L","P"))
erosion_dat$Seeding <- as.factor(erosion_dat$Seeding)
erosion_dat$BSNE_gM2day <- as.numeric(erosion_dat$BSNE_gM2day)

BSNE_dat <- erosion_dat %>%
  filter(BSNE_gM2day != "#N/A") #create new wind erosion data set with the the NA values removed

str(BSNE_dat)

SMBSNE_dat <- subset(BSNE_dat, Site=="SM") #Create new data sets by site so I run sepearate models by site
WMBSNE_dat <- subset(BSNE_dat, Site=="WM")

#Remove NAs from Chla data for the Silt_erosion analysis
erosion_dat <- erosion_dat %>%
  filter(Chla != "NA")

###Explore a simple lmer fit with random variables for Wind Erosion 
lmerfit.Wind <- glmer((BSNE_gM2day+.01) ~ Treatment_Year + VegCover + Median_SoilStability + Soil + Chla + (1|Transect) + (1|Treatment) + (1|Seeding) + (1|Site), family=Gamma(link="log"), data=BSNE_dat)
summary(lmerfit.Wind)
plot(lmerfit.Wind) #This looks terrible
qqnorm(residuals(lmerfit.Wind)) #This doesn't look bad
#Error that model failed to converge. 

#Check singularity
tt <- getME(lmerfit.Wind,"theta")
ll <- getME(lmerfit.Wind,"lower")
min(tt[ll==0]) #Don't want this to be too close to zero; it is 0.0005591668. This could be an issue, but it could be fine.

#Double check gradient calculations
derivs1 <- lmerfit.Wind@optinfo$derivs
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1)) #This should be <0.001, and it is: 5.192835e-05

#Retry convergence and bump up the number of iterations?
ss <- getME(lmerfit.Wind,c("theta","fixef"))
m2 <- update(lmerfit.Wind,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4))) #this did not work at all
#Trying a different optimizer...bobyqa
m3 <- update(lmerfit.Wind,start=ss,control=glmerControl(optimizer="Nelder_Mead",
                                                 optCtrl=list(maxfun=2e5))) #This did not converge; degenerate Hessian with 1 neg eigen values
m4 <- update(lmerfit.Wind,start=ss,control=glmerControl(optimizer="bobyqa",
                                                        optCtrl=list(maxfun=2e5)))
#I am stuck. Not sure how to fix this problem


SMlmerfit.Wind <- glmer((BSNE_gM2day+.01) ~ Treatment_Year + VegCover + Median_SoilStability + Soil + Chla + (1|Transect) + (1|Treatment) + (1|Seeding), family=Gamma(link="log"), data=SMBSNE_dat)
summary(SMlmerfit.Wind)
plot(SMlmerfit.Wind) #This looks okay I think
qqnorm(residuals(SMlmerfit.Wind)) #This does not look good

WMlmerfit.Wind <- glmer((BSNE_gM2day+.01) ~ Treatment_Year + VegCover + Median_SoilStability + Soil + Chla + (1|Transect) + (1|Treatment) + (1|Seeding), family=Gamma(link="log"), data=WMBSNE_dat)
summary(WMlmerfit.Wind) 
plot(WMlmerfit.Wind) #This looks bad
qqnorm(residuals(WMlmerfit.Wind)) #not bad!

#Neither of these model converge. I don't really know what that means. I still get parameter estimates for them.

###GLMMs didn't converge in glmer() so I'm moving on to a Bayesian approach...
library(tidyr)
library(rstan)     #for extract()
library(rstanarm)  #Bayesian multilevel: stan_lmer(), stan_glmer() etc
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
theme_set(theme_grey()) #rstanarm overrides default ggplot theme: set it back
source("hpdi.R") #For calculating credible intervals  

#Fit a Bayesian model for Wind Erosion
bysfitWind <- stan_glmer((BSNE_gM2day+.01) ~ Treatment_Year + VegCover + Median_SoilStability + Soil + Chla + (1|Transect) + (1|Treatment) + (1|Seeding) + (1|Site),
                      family=Gamma(link="log"), data=BSNE_dat)
summary(bysfitWind,digits=4)

vcov(bysfitWind,correlation=TRUE) #Correlation matrix
bysfitWindsamples <- extract(bysfitWind$stanfit) #extract samples
names(bysfitWindsamples)
hist(bysfitWindsamples$alpha[,1],breaks=75,col="lightblue") #\beta_0 = intercept = alpha
hist(bysfitWindsamples$beta[,1],breaks=75,col="lightblue") #\beta_1
hist(bysfitWindsamples$beta[,3],breaks=75,col="lightblue") #\beta_3

#Fit another Bayesian Wind Erosion Model with nested random effects

bysfitWind.2 <- stan_glmer(log(BSNE_gM2day+.01) ~ Treatment_Year + VegCover + Median_SoilStability + Soil + Chla + (1|Site/Seeding/Treatment/Transect),
                         data=BSNE_dat)
summary(bysfitWind.2)
vcov(bysfitWind.2,correlation=TRUE) #Nothing is too correlated. 
samples.BW2 <- extract(bysfitWind.2$stanfit)
names(samples.BW2)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

samplesdf.BW2 <- data.frame(samples.BW2$alpha,samples.BW2$beta)
samplesdf.BW2
names(samplesdf.BW2) <- c("alpha",paste(names(samples.BW2[2]),1:5,sep="_"))

hist(samples.BW2$alpha[,1],breaks=75,col="lightblue") #\beta_0 = intercept = alpha
hist(samples.BW2$beta[,1],breaks=75,col="lightblue") #\beta_1
hist(samples.BW2$beta[,2],breaks=75,col="lightblue") #\beta_2
hist(samples.BW2$beta[,3],breaks=75,col="lightblue") #\beta_3
hist(samples.BW2$beta[,4],breaks=75,col="lightblue") #\beta_4
hist(samples.BW2$beta[,5],breaks=75,col="lightblue") #\beta_5
#All of these histograms look pretty normal

launch_shinystan(bysfitWind.2)
#That looks pretty good too!!

#Now, run the Wind erosion models separated by site
bysfitWind.SM <- stan_glmer(log(BSNE_gM2day+.01) ~ Treatment_Year + VegCover + Median_SoilStability + Soil + Chla + (1|Seeding/Treatment/Transect),
                           data=SMBSNE_dat)
bysfitWind.WM <- stan_glmer(log(BSNE_gM2day+.01) ~ Treatment_Year + VegCover + Median_SoilStability + Soil + Chla + (1|Seeding/Treatment/Transect),
                            data=WMBSNE_dat)

vcov(bysfitWind.SM,correlation=TRUE)
samples.BSM <- extract(bysfitWind.SM$stanfit)
samplesdf.BSM <- data.frame(samples.BSM$alpha,samples.BSM$beta)
names(samplesdf.BSM) <- c("alpha",paste(names(samples.BSM[2]),1:5,sep="_"))
samplesdf.BSM

