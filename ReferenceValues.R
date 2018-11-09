reference_dat <- read.csv("ReferenceValues.csv", as.is=TRUE, na.strings = "n/a") #read in data

head(reference_dat)

meansVegCover <- aggregate(reference_dat[, 8:14], list(reference_dat$Site), mean)
SMmeanSoilCover <- aggregate(SM[, 15:27], list(SM$Site), mean)


SM <- subset(reference_dat, Site == "SM")
str(SM)

colMeans(SM[sapply(SM, is.character)])

calcSE <- function(x){sd(x)/sqrt(length(x))}
VegCoverSE <- aggregate(reference_dat[, 8:14], list(reference_dat$Site), calcSE)
VegCoverSE
