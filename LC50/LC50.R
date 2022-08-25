rm(list =ls())


install.packages('drc')


library(ggplot2)
library(dplyr)
###### Disables scientific notation #####
options(scipen = 999)
library(drc)
Yes

###Data
D <- read.csv("LC50 R.csv")

Imidacloprid5 <- filter(D, D$Treatment=='imidacloprid', D$Temperature.ºC == 5)
print(Imidacloprid5)

###################################Fipronil 5ºC 
F5 <- read.csv("Fipronil 5ºC.csv.csv")

#Calculate LC50 F5
LC50_F5 <- drm(Dead/Total~Concentration,weights=Total,fct=LL.2(),type="binomial",data=F5)
summary(LC50_F5)
plot(LC50_F5, xlab="Concentration (µg/L)", ylab="mortality (%)")
ED(LC50_F5,c(10,50,90))


###################################Imidacloprid 5ªC
I5 <- read.csv("Imidacloprid 5ºC.csv")

#Calculate LC50 I5
LC50_I5 <- drm(Dead/Total~Concentration,weights=Total,fct=LL.2(),type="binomial",data=I5)
summary(LC50_I5)
plot(LC50_I5)
ED(LC50_I5,c(10,50,90))


###############################Imidacloprid 15
I15 <- read.csv("Imidacloprid 15ºC.csv")

#Calculate LC50 I15
LC50_I15 <- drm(Dead/Total~Concentration,weights=Total,fct=LL.2(),type="binomial",data=I15)
summary(LC50_I15)
plot(LC50_I15, xlab="Concentration (µg/L)", ylab="mortality (%)")
ED(LC50_I15,c(10,50,90))


##############Fluralaner 5ºC
Fl5 <- read.csv("Fluralaner 5ºC.csv")

#Calculate LC50 Fl5
LC50_Fl5 <- drm(Dead/Total~Concentration,weights=Total,fct=LL.2(),type="binomial",data=Fl5)
summary(LC50_Fl5)
plot(LC50_Fl5,xlab="Concentration (µg/L)", ylab="mortality (%)")
ED(LC50_Fl5,c(10,50,90))

###################Selamectin 5ºC
S5  <- read.csv("Selamectin 5ºC.csv")


#Calculate LC50 S5
LC50_S5 <- drm(Dead/Total~Concentration,weights=Total,fct=LL.2(),type="binomial",data=S5)
summary(LC50_S5)
plot(LC50_S5, xlab="Concentration (µg/L)", ylab="mortality (%)")
ED(LC50_S5,c(10,50,90))








