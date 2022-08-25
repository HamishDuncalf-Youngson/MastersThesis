rm(list=ls())


d <- read.csv("GammarusMobilityTrials.csv")
str(d)


library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)

# Function gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



#Barchart for all Distance by temperature 
distancetemperature<- summarySE(d,measurevar= "Distance", groupvars=c("Temperature"),na.rm = TRUE)
distancetemperature2 <- distancetemperature
distancetemperature2$Temperature <- factor(distancetemperature2$Temperature)
ggplot(distancetemperature2, aes(x= Temperature, y=Distance, fill=Temperature)) +   geom_bar(position = "dodge", stat = "identity") + 
  geom_errorbar(aes(ymin= Distance-se, ymax= Distance+se),  width=.22, position=position_dodge(.9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs (x= "Temperature (ºC)",  y="Distance (cm)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))




#Barchart for all Distance by chemical
distancechemical<- summarySE(d,measurevar= "Distance", groupvars=c("Chemical"),na.rm = TRUE)
distancechemical2 <- distancechemical
distancechemical2$Chemical <- factor(distancechemical2$Chemical)
ggplot(distancechemical2, aes(x= Chemical, y=Distance, fill=Chemical)) +   geom_bar(position = "dodge", stat = "identity") + 
  geom_errorbar(aes(ymin= Distance-se, ymax= Distance+se),  width=.22, position=position_dodge(.9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs (x= "Chemical Treatment",  y="Distance (cm)")  + scale_fill_lancet()



#Scatterplot of Distance by body length 
plot(d$Distance~d$Body.Length, ylab="Distance (cm)",xlab="Body Length (cm)")
abline(m)

#linear model
m <- lm(d$Distance~d$Body.Length)
summary(m)
abline(m)

##Histogram of body lengths
hist(d$Body.Length)


#SE for the whole Data
seDistance<-sqrt(var(d$Distance)/length(d$Distance)) 
seDistance
mean((d$Distance))

#SE for 10ºC
d10<-subset(d, d$Temperature==10)
SEd10<-sqrt(var(d10$Distance)/length(d10$Distance))
SEd10
mean(d10$Distance)

#SE for 15ºC
d15<-subset(d, d$Temperature==15)
SEd15<-sqrt(var(d15$Distance)/length(d15$Distance))
SEd15
mean(d15$Distance)

#SE for control
dCntrl<-subset(d, d$Chemical== 'No Chemical')
SEdCntrl<-sqrt(var(dCntrl$Distance)/length(dCntrl$Distance))
SEdCntrl
mean(dCntrl$Distance)

#SE for Imidacloprid
dImida<-subset(d, d$Chemical== 'Imidacloprid')
SEdImida<-sqrt(var(dImida$Distance)/length(dImida$Distance))
SEdImida
mean(dImida$Distance)


hist(d$Distance,breaks=213)
hist(d10$Distance)
hist(d15$Distance)
hist(dCntrl$Distance)
hist(dImida$Distance)

hist(log(d$Distance+1),breaks=20,col= '#35DAC1')
hist(log(d10$Distance+1),breaks=20,col= '#35DAC1')
hist(log(d15$Distance+1),breaks=20,col= '#35DAC1')
hist(log(dCntrl$Distance+1),breaks=20,col= '#35DAC1')
hist(log(dImida$Distance+1),breaks=20, col= '#35DAC1')

par(mfrow=c(1,1))


logd10<- log(d10$Distance+1)
logd <- log(d$Distance+1)
logd2 <- log(d$Distance)
logv <- log(d$Variance)

logd
logd2


boxplot(logd)

SEtest<-sqrt(var(logd10)/length(logd10))
SEtest
mean(logd10)


test2 <- rbind(d,logd)
View(test2)


library(dplyr)

#Workss for combining log data with original

f <- mutate(d,Log_Variance=logv)
d <- mutate(f,Log_Distance=logd2)

##Standard error for total distance

#SE for 5ºC Control
d5control<-filter(d, d$Chemical== 'No Chemical', d$Temperature==5)
SEd5control<-sqrt(var(d5control$Distance)/length(d5control$Distance))
SEd5control
mean(d5control$Distance)

#SE for 10ºC Control
d10Control<-filter(d10, d10$Chemical== 'No Chemical')
SEd10Control<-sqrt(var(d10Control$Distance)/length(d10Control$Distance))
SEd10Control
mean(d10Control$Distance)

#SE for 15ºC Control 
d15Control<- filter(d15, d15$Chemical== 'No Chemical')
SEd15Control<-sqrt(var(d15Control$Distance)/length(d15Control$Distance))
SEd15Control
mean(d15Control$Distance)

#SE for 5ºC Imidacloprid
d5imida<-filter(d, d$Chemical== 'Imidacloprid', d$Temperature==5)
SEd5imida<-sqrt(var(d5imida$Distance)/length(d5imida$Distance))
SEd5imida
mean(d5imida$Distance)

#SE for 10ºC Imidacloprid 
d10Imida<-filter(d10, d10$Chemical== 'Imidacloprid')
SEd10Imida<-sqrt(var(d10Imida$Distance)/length(d10Imida$Distance))
SEd10Imida
mean(d10Imida$Distance)

#SE for 15ºC Imidacloprid
d15Imida<- filter(d15, d15$Chemical== 'Imidacloprid')
SEd15Imida<-sqrt(var(d15Imida$Distance)/length(d15Imida$Distance))
SEd15Imida
mean(d15Imida$Distance)

#SE for 5ºC Pesticide Mix 1 5ºC ###
f5p1<-filter(c, c$Trial== 'P Mix 2 15ºC')
SEFP215<-sqrt(var(FP215$Distance.cm)/length(FP215$Distance.cm))
SEFP215
mean(FP215$Distance.cm)


##############################################
#Subsetting for log ((post event?) turns out I didnt need to do this, see below))


f10<-subset(f, f$Temperature==10)
f15<-subset(f, f$Temperature==15)
fCntrl<-subset(f, f$Chemical== 'No Chemical')
fImida<-subset(f, f$Chemical== 'Imidacloprid')


#################################################
##Standard error of log total distance

#SE for log 5ºC Control 
f5control<- filter(d, d$Chemical== 'No Chemical',d$Temperature==5)
SEf5Control<-sqrt(var(f5control$Log_Distance)/length(f5control$Log_Distance))
SEf5Control
mean(f5control$Log_Distance)

#SE for log 10ºC Control
f10Control<-filter(f10, f10$Chemical== 'No Chemical')
SEf10Control<-sqrt(var(f10Control$Log_Distance)/length(f10Control$Log_Distance))
SEf10Control
mean(f10Control$Log_Distance)

#SE for log 15ºC Control 
f15Control<- filter(f15, f15$Chemical== 'No Chemical')
SEf15Control<-sqrt(var(f15Control$Log_Distance)/length(f15Control$Log_Distance))
SEf15Control
mean(f15Control$Log_Distance)

#SE for log 5ºC Imidacloprid 
f5Imida<- filter(d, d$Chemical== 'Imidacloprid',d$Temperature==5)
SEf5Imida<-sqrt(var(f5Imida$Log_Distance)/length(f5Imida$Log_Distance))
SEf5Imida
mean(f5Imida$Log_Distance)

#SE for log 10ºC Imidacloprid 
f10Imida<-filter(f10, f10$Chemical== 'Imidacloprid')
SEf10Imida<-sqrt(var(f10Imida$Log_Distance)/length(f10Imida$Log_Distance))
SEf10Imida
mean(f10Imida$Log_Distance)

#SE for log 15ºC Imidacloprid
f15Imida<- filter(f15, f15$Chemical== 'Imidacloprid')
SEf15Imida<-sqrt(var(f15Imida$Log_Distance)/length(f15Imida$Log_Distance))
SEf15Imida
mean(f15Imida$Log_Distance)


##############*Remember Log Boxplot*################################


boxplot(f$Log_Distance~d$Temperature+d$Chemical, 
        col = c("orchid2", "orchid2","#1AD6AE","#1AD6AE"), ylab="Log Distance (cm)", xlab= "")


boxplot(f$Distance~f$Temperature+f$Chemical, 
        col = c("orchid2", "orchid2","#1AD6AE","#1AD6AE"), ylab="Distance (cm)", xlab= "")

v <- read.csv("Gammarus_Steps.csv")
v1 <- read.csv("Gammarus_steps2.csv")

View(v1)


#Trying to do some variance stuff
var(v$Gam.1.Distance)
mean(v$Gam.1.Distance)
sum(v$Gam.1.Distance)


hmm <- list(13,48,58)



attach(v)
hello <- (c((var(Gam.1.Distance)),(var(Gam.2.Distance)),
           (var(Gam.3.Distance)),(var(Gam.4.Distance)),(var(Gam.5.Distance)),(var(Gam.6.Distance)),(var(Gam.7.Distance)),
           (var(Gam.8.Distance)),(var(Gam.9.Distance)),(var(Gam.10.Distance)),(var(Gam.11.Distance)),(var(Gam.12.Distance))
           ,(var(Gam.13.Distance)),(var(Gam.14.Distance)),(var(Gam.15.Distance)),(var(Gam.16.Distance)),(var(Gam.17.Distance))
           ,(var(Gam.18.Distance)),(var(Gam.19.Distance)),(var(Gam.20.Distance)),(var(Gam.21.Distance)),(var(Gam.22.Distance))
           ,(var(Gam.23.Distance)),(var(Gam.24.Distance))))

attach(v)           
Goodbye <- (c((mean(Gam.1.Distance)),(mean(Gam.2.Distance)),
            (mean(Gam.3.Distance)),(mean(Gam.4.Distance)),(mean(Gam.5.Distance)),(mean(Gam.6.Distance)),(mean(Gam.7.Distance)),
            (mean(Gam.8.Distance)),(mean(Gam.9.Distance)),(mean(Gam.10.Distance)),(mean(Gam.11.Distance)),(mean(Gam.12.Distance))
            ,(mean(Gam.13.Distance)),(mean(Gam.14.Distance)),(mean(Gam.15.Distance)),(mean(Gam.16.Distance)),(mean(Gam.17.Distance))
            ,(mean(Gam.18.Distance)),(mean(Gam.19.Distance)),(mean(Gam.20.Distance)),(mean(Gam.21.Distance)),(mean(Gam.22.Distance))
            ,(mean(Gam.23.Distance)),(mean(Gam.24.Distance))))

str(hello)
str(Goodbye)
#Wahey
maybe <- mutate(Bug,Variance= hello)
maybenot <- mutate(maybe,Mean= Goodbye)
View(maybenot)



boxplot(maybenot$Variance~maybenot$Temperature+maybenot$Chemical, 
        col = c("orchid2", "orchid2","#1AD6AE","#1AD6AE"), ylab="Variance", xlab= "")

boxplot(maybenot$Mean~maybenot$Temperature+maybenot$Chemical, 
        col = c("orchid2", "orchid2","#1AD6AE","#1AD6AE"), ylab="Mean Step (cm)", xlab= "")

### Standard Error of the variance#######

#SE for Variance 5ºC Control 
fv5Control<-filter(d, d$Chemical== 'No Chemical', d$Temperature==5)
SEfv5Control<-sqrt(var(fv5Control$Variance)/length(fv5Control$Variance))
SEfv5Control
mean(fv5Control$Variance)

#SE for Variance 10ºC Control
fv10Control<-filter(maybenot, maybenot$Chemical== 'No Chemical', maybenot$Temperature==10)
SEfv10Control<-sqrt(var(fv10Control$Variance)/length(fv10Control$Variance))
SEfv10Control
mean(fv10Control$Variance)

#SE for log 15ºC Control 
fv15Control<- filter(maybenot, maybenot$Chemical== 'No Chemical',maybenot$Temperature==15)
SEfv15Control<-sqrt(var(fv15Control$Variance)/length(fv15Control$Variance))
SEfv15Control
mean(fv15Control$Log_Distance)

#SE for Variance 5ºC Imidacloprid 
fv5Imida<-filter(d, d$Chemical== 'Imidacloprid', d$Temperature==5)
SEfv5Imida<-sqrt(var(fv5Imida$Variance)/length(fv5Imida$Variance))
SEfv5Imida
mean(fv5Imida$Variance)

#SE for log 10ºC Imidacloprid 
fv10Imida<-filter(maybenot, maybenot$Chemical== 'Imidacloprid', maybenot$Temperature==10)
SEfv10Imida<-sqrt(var(fv10Imida$Variance)/length(fv10Imida$Variance))
SEfv10Imida
mean(fv10Imida$Variance)

#SE for log 15ºC Imidacloprid
fv15Imida<- filter(maybenot, maybenot$Chemical== 'Imidacloprid', maybenot$Temperature==15)
SEvf15Imida<-sqrt(var(fv15Imida$Variance)/length(fv15Imida$Variance))
SEvf15Imida
mean(fv15Imida$Variance)


##Boxplot##

boxplot(maybeyes$logvariance~maybeyes$Temperature+maybeyes$Chemical, 
        col = c("orchid2", "orchid2","#1AD6AE","#1AD6AE"), ylab="Log Variance", xlab= "")

write.csv(maybenot,"variance.csv")

#Maybe log variance next?!

vlog <- log(maybenot$Variance)
vlog


maybeyes <- mutate(maybenot,logvariance=vlog)


write.csv(maybeyes,"variancelog.csv")


##Standard error of log Variance##

##SE for log Variance 5ºC Control##
fvl5Control<-filter(d, d$Chemical== 'No Chemical', d$Temperature==5)
SEfvl5Control<-sqrt(var(fvl5Control$Log_Variance)/length(fvl5Control$Log_Variance))
SEfvl5Control
mean(fvl5Control$Log_Variance)


##SE for log Variance 10ºC Control##
fvl10Control<-filter(maybeyes, maybeyes$Chemical== 'No Chemical', maybeyes$Temperature==10)
SEfvl10Control<-sqrt(var(fvl10Control$logvariance)/length(fvl10Control$logvariance))
SEfvl10Control
mean(fvl10Control$logvariance)

#SE for log 15ºC Control 
fvl15Control<- filter(maybeyes, maybeyes$Chemical== 'No Chemical',maybeyes$Temperature==15)
SEfvl15Control<-sqrt(var(fvl15Control$logvariance)/length(fvl15Control$logvariance))
SEfvl15Control
mean(fvl15Control$logvariance)

#SE for log 5ºC Imidacloprid
fvl5Imida<-filter(d, d$Chemical== 'Imidacloprid', d$Temperature==5)
SEfvl5Imida<-sqrt(var(fvl5Imida$Log_Variance)/length(fvl5Imida$Log_Variance))
SEfvl5Imida
mean(fvl5Imida$Log_Variance)

#SE for log 10ºC Imidacloprid 
fvl10Imida<-filter(maybeyes, maybeyes$Chemical== 'Imidacloprid', maybeyes$Temperature==10)
SEfvl10Imida<-sqrt(var(fvl10Imida$logvariance)/length(fvl10Imida$logvariance))
SEfvl10Imida
mean(fvl10Imida$logvariance)

#SE for log 15ºC Imidacloprid
fvl15Imida<- filter(maybeyes, maybeyes$Chemical== 'Imidacloprid', maybeyes$Temperature==15)
SEfvl15Imida<-sqrt(var(fvl15Imida$logvariance)/length(fvl15Imida$logvariance))
SEfvl15Imida
mean(fvl15Imida$logvariance)




boxplot(d$Distance~d$Temperature+d$Chemical, col = c("orchid2", "#1AD6AE","orchid2","#1AD6AE"), ylab="Distance (cm)", xlab= "")

install.packages("ggplot2")
library(ggplot2)



ggplot(c, aes(x = Trial y = Distance, fill = Trial )) + 
        geom_boxplot() 


c <- read.csv("Mobilitymaster.csv")

get.se <- function(y) {
    se <- sd(y)/sqrt(length(y))
    mu <- mean(y)
    c(ymin=mu-se, ymax=mu+se)
}


###Barchart for total distance Imidacloprid###
###and whichever other ones I want to compare to control.... 

c$cont<-c(rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6))
c$cont<-factor(c$cont, levels=c("Control", "Exp"))

trialdistance<- summarySE(c,measurevar= "Distance.cm", groupvars=c("Trial"),na.rm = TRUE)
trialdistance2 <- trialdistance
trialdistance2$Trial <- factor(trialdistance2$Trial)

trialdistance2$chemical <- c(rep("Control", 3),rep("Fipronil",3),rep("Fluralaner",3), rep("Imidacloprid",3), rep("P Mix 1", 3), rep("P Mix 2", 3), rep("Selamectin",3))

ggplot(trialdistance2, aes(x=Trial, y=Distance.cm, fill=chemical)) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(title = "Imidacloprid 100 (µg/L)", x= "",  y="Distance (cm)") +
  geom_errorbar(aes(ymin= Distance.cm-se, ymax= Distance.cm+se),  width=.22, position=position_dodge(.9)) +
  scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC"))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none")
    
write.csv(d,"Allthings.csv")





###Faceted all chems at 100µg/L

trialdistance2$Trial2 <- factor(trialdistance2$Trial, levels = c("Control 5ºC","Control 10ºC","Control 15ºC","Fipronil 5ºC ","Fipronil 10ºC", "Fipronil 15ºC ", 
                           "Fluralaner 5ºC ", "Fluralaner 10ºC ", "Fluralaner 15ºC ", "Imidacloprid 5ºC", "Imidacloprid 10ºC", "Imidacloprid 15ºC", "P Mix 1 5ºC ", "P Mix 1 10ºC", "P Mix 1 15ºC", 
                           "P Mix 2 5ºC", "P Mix 2 10ºC ", "P Mix 2 15ºC", "Selamectin 5ºC", "Selamectin 10ºC", "Selamectin 15ºC "))

ggplot(trialdistance2, aes(x=Trial2, y=Distance.cm, fill=chemical)) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(title = "Imidacloprid 100 (µg/L)", x= "",  y="Distance (cm)") +
  geom_errorbar(aes(ymin= Distance.cm-se, ymax= Distance.cm+se),  width=.22, position=position_dodge(.9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none") + 
  facet_wrap(~chemical, scales='free_x') + 
  scale_fill_jama()
                                     



 
###Barchart for log total distance###

c <- read.csv("Mobilitymaster.csv")


c$cont<-c(rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6))
c$cont<-factor(c$cont, levels=c("Control", "Exp"))
ggplot(c, aes(x=Trial, y= Log.Distance, fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "turquoise"))+
    theme(legend.position = "none")



##Barchart for Variance##
a <- read.csv("Variancesheet2.csv")


a$cont<-c(rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6))
a$cont<-factor(a$cont, levels=c("Control", "Exp"))
ggplot(a, aes(x=Trial, y=Variance, fill=cont)) +
        stat_summary(fun=mean, geom="bar")+
        stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
        scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC"))+
        scale_fill_manual(values=c("mediumorchid2", "turquoise"))+
        theme(legend.position = "none")


##Barchart for log Variance##
write.csv(a,"Variancesheet2.csv")
a <- read.csv("Variancesheet2.csv")


a$cont<-c(rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6))
a$cont<-factor(a$cont, levels=c("Control", "Exp"))
ggplot(a, aes(x=Trial, y=Log.Variance, fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "turquoise"))+
    theme(legend.position = "none")


##Barchart for total Area##
ggplot(c, aes(x=Trial, y=Total.Area, fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "turquoise"))+
    theme(legend.position = "none")



#Total Area Log
Tlog <- log(c$Total.Area)
Tlog


ttest <- mutate(c,TotalAreaLog=Tlog)
c <- ttest

write.csv(ttest,"TAlog.csv")

##Barchart for total Area log##
ggplot(c, aes(x=Trial, y=TotalAreaLog, fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "turquoise"))+
    theme(legend.position = "none")


###Barchart for 24hr###

t24 <- read.csv("24hrs/24Master.csv")


t24$cont<-c(rep("Control", 6), rep("Exp", 6),
            rep("Control", 6), rep("Exp", 6),
            rep("Control", 6), rep("Exp", 6))
t24$cont<-factor(a$cont, levels=c("Control", "Exp"))
ggplot(t24, aes(x=Trial, y=Distance.cm..24hour., fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "turquoise"))+
    theme(legend.position = "none")



###Barchart for 24hr Variance###

t24 <- read.csv("24hrs/24Master.csv")


t24$cont<-c(rep("Control", 6), rep("Exp", 6),
            rep("Control", 6), rep("Exp", 6),
            rep("Control", 6), rep("Exp", 6))
t24$cont<-factor(a$cont, levels=c("Control", "Exp"))
ggplot(t24, aes(x=Trial, y= Variance, fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "turquoise"))+
    theme(legend.position = "none")



### Log for 24hr total distacnce##
t24log <- log(t24$Distance.cm..24hour.+1)
t24log

t24test <- mutate(t24, 'Distance.cm.24hourlog'= t24log)
t24 <- t24test

write.csv(t24,"24hourmaster2.csv")



### Log for 24hr Variance##
t24logvar <- log(t24$Variance+1)
t24logvar

t24test <- mutate(t24, variancelog = t24logvar)
t24 <- t24test

write.csv(t24,"24hourmaster3.csv")

##Barchart for 24hourlog distance ##

ggplot(t24, aes(x=Trial, y=Distance.cm.24hourlog, fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "turquoise"))+
    theme(legend.position = "none") +labs (x= "",y="Distance log (x+1)")

##Barchart for 24hourlog Variance##

ggplot(t24, aes(x=Trial, y=variancelog, fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "turquoise"))+
    theme(legend.position = "none") +labs (x= "",y="Variance Log (x+1)")

#SE for 24hour 5ºC Control 
f245c<-filter(t24, t24$Trial== 'Control 5ºC')
SEf245c<-sqrt(var(f245c$Distance.cm..24hour.)/length(f245c$Distance.cm..24hour.))
SEf245c
mean(f245c$Distance.cm..24hour.)

#SE for 24hour 10ºC Control
f2410c<-filter(t24, t24$Trial== 'Control 10ºC')
SEf2410c<-sqrt(var(f2410c$Distance.cm..24hour.)/length(f2410c$Distance.cm..24hour.))
SEf2410c
mean(f2410c$Distance.cm..24hour.)

#SE for 24hour 15ºC Control 
f2415c<-filter(t24, t24$Trial== 'Control 15ºC')
SEf2415c<-sqrt(var(f2415c$Distance.cm..24hour.)/length(f2415c$Distance.cm..24hour.))
SEf2415c
mean(f2415c$Distance.cm..24hour.)

#SE for 24hour 5ºC Imidacloprid 
f245i<-filter(t24, t24$Trial== 'Imidacloprid 5ºC')
SEf245i<-sqrt(var(f245i$Distance.cm..24hour.)/length(f245i$Distance.cm..24hour.))
SEf245i
mean(f245i$Distance.cm..24hour.)

#SE for 24hour 10ºC Imidacloprid 
f2410i<-filter(t24, t24$Trial== 'Imidacloprid 10ºC')
SEf2410i<-sqrt(var(f2410i$Distance.cm..24hour.)/length(f2410i$Distance.cm..24hour.))
SEf2410i
mean(f2410i$Distance.cm..24hour.)

#SE for 24hour 15ºC Imidacloprid 
f2415i<-filter(t24, t24$Trial== 'Imidacloprid 15ºC')
SEf2415i<-sqrt(var(f2415i$Distance.cm..24hour.)/length(f2415i$Distance.cm..24hour.))
SEf2415i
mean(f2415i$Distance.cm..24hour.)

###Barchart for pesticide mix 1###

x <- read.csv("MobilityMasterPesticide1.csv")
c <- read.csv("Mobilitymaster.csv")

x$cont<-c(rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6))
x$cont<-factor(x$cont, levels=c("Control", "Exp"))
ggplot(x, aes(x=Trial, y=Distance.cm, fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC", "P Mix 1 5ºC ", "P Mix 1 10ºC", "P Mix 1 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "seagreen3"))+
    theme(legend.position = "none")


###Barchart for pesticide mix 2###


y <- read.csv("MobilityMasterPesticide2.csv")

y$cont<-c(rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6))
y$cont<-factor(x$cont, levels=c("Control", "Exp"))
ggplot(y, aes(x=Trial, y=Distance.cm, fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC", "P Mix 2 5ºC", "P Mix 2 10ºC ", "P Mix 2 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "goldenrod2"))+
    theme(legend.position = "none")



##Barchart for all ###

c$cont<-c(rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Exp1", 12), rep("Exp1", 12),
          rep("Exp1", 12), rep("Exp2", 12),
          rep("Exp2", 12), rep("Exp2", 12))
c$cont<-factor(c$cont, levels=c("Control", "Exp", "Exp1", "Exp2"))
ggplot(c, aes(x=Trial, y=Distance.cm, fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC","P Mix 1 5ºC ", "P Mix 1 10ºC", "P Mix 1 15ºC", "P Mix 2 5ºC", "P Mix 2 10ºC ", "P Mix 2 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "turquoise","goldenrod2", "seagreen3"))+
    theme(legend.position = "none", panel.background = element_blank())


### Log for all mix ###

allLog1 <- log(c$Distance.cm+1)
allLog1

alllogtest <- mutate(c, 'TotalLogDistance1'= allLog1)
c <- alllogtest

##Barchart for all log ##

c$cont<-c(rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Exp1", 12), rep("Exp1", 12),
          rep("Exp1", 12), rep("Exp2", 12),
          rep("Exp2", 12), rep("Exp2", 12))
c$cont<-factor(c$cont, levels=c("Control", "Exp", "Exp1", "Exp2"))
ggplot(c, aes(x=Trial, y=TotalLogDistance1 , fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC","P Mix 1 5ºC ", "P Mix 1 10ºC", "P Mix 1 15ºC", "P Mix 2 5ºC", "P Mix 2 10ºC ", "P Mix 2 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "turquoise","goldenrod2", "seagreen3"))+
    theme(legend.position = "none", panel.background = element_blank())  +labs (x= "",y="Distance.cm Log (x+1)")


write.csv(c,"MobilityMaster3.csv")

##Barchart for all variance###

c$cont<-c(rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Exp1", 12), rep("Exp1", 12),
          rep("Exp1", 12), rep("Exp2", 12),
          rep("Exp2", 12), rep("Exp2", 12))
c$cont<-factor(c$cont, levels=c("Control", "Exp", "Exp1", "Exp2"))
ggplot(c, aes(x=Trial, y=Variance, fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC","P Mix 1 5ºC ", "P Mix 1 10ºC", "P Mix 1 15ºC", "P Mix 2 5ºC", "P Mix 2 10ºC ", "P Mix 2 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "turquoise","goldenrod2", "seagreen3"))+
    theme(legend.position = "none", panel.background = element_blank())  +labs (x= "",y="Variance")

### Log for all variance ###

allLog1var <- log(c$Variance+1)
allLog1var

alllogtestvar <- mutate(c, 'Variancelog'= allLog1var)
c <- alllogtestvar

##Barchart for all log variance ##

c$cont<-c(rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Exp1", 12), rep("Exp1", 12),
          rep("Exp1", 12), rep("Exp2", 12),
          rep("Exp2", 12), rep("Exp2", 12))
c$cont<-factor(c$cont, levels=c("Control", "Exp", "Exp1", "Exp2"))
ggplot(c, aes(x=Trial, y= Variancelog , fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC","P Mix 1 5ºC ", "P Mix 1 10ºC", "P Mix 1 15ºC", "P Mix 2 5ºC", "P Mix 2 10ºC ", "P Mix 2 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "turquoise","goldenrod2", "seagreen3"))+
    theme(legend.position = "none", panel.background = element_blank())  +labs (x= "",y="Variance Log (x+1)")


###SE for Pesticide Mix 1 5ºC ###
FP15<-filter(c, c$Trial== 'P Mix 1 5ºC ')
SEFP15<-sqrt(var(FP15$Distance.cm)/length(FP15$Distance.cm))
SEFP15
mean(FP15$Distance.cm)

###SE for Pesticide Mix 1 10ºC ###
FP110<-filter(c, c$Trial== 'P Mix 1 10ºC')
SEFP110<-sqrt(var(FP110$Distance.cm)/length(FP110$Distance.cm))
SEFP110
mean(FP110$Distance.cm)

###SE for Pesticide Mix 1 15ºC ###
FP115<-filter(c, c$Trial== 'P Mix 1 15ºC')
SEFP115<-sqrt(var(FP115$Distance.cm)/length(FP115$Distance.cm))
SEFP115
mean(FP115$Distance.cm)

###SE for Pesticide Mix 2 5ºC ###
FP25<-filter(c, c$Trial== 'P Mix 2 5ºC')
SEFP25<-sqrt(var(FP25$Distance.cm)/length(FP25$Distance.cm))
SEFP25
mean(FP25$Distance.cm)

###SE for Pesticide Mix 2 10ºC ###
FP210<-filter(c, c$Trial== 'P Mix 2 10ºC ')
SEFP210<-sqrt(var(FP210$Distance.cm)/length(FP210$Distance.cm))
SEFP210
mean(FP210$Distance.cm)

###SE for Pesticide Mix 2 15ºC ###
FP215<-filter(c, c$Trial== 'P Mix 2 15ºC')
SEFP215<-sqrt(var(FP215$Distance.cm)/length(FP215$Distance.cm))
SEFP215
mean(FP215$Distance.cm)




#### SE's for logged Pesticide Mixes ###

###SE for log Pesticide Mix 1 5ºC ###
FP15L<-filter(c, c$Trial== 'P Mix 1 5ºC ')
SEFP15L<-sqrt(var(FP15L$TotalLogDistance1)/length(FP15$TotalLogDistance1))
SEFP15L
mean(FP15L$TotalLogDistance1)


###SE for log Pesticide Mix 1 10ºC ###
FP110L<-filter(c, c$Trial== 'P Mix 1 10ºC')
SEFP110L<-sqrt(var(FP110L$TotalLogDistance1)/length(FP110L$TotalLogDistance1))
SEFP110L
mean(FP110L$TotalLogDistance1)

###SE for log Pesticide Mix 1 15ºC ###
FP115L<-filter(c, c$Trial== 'P Mix 1 15ºC')
SEFP115L<-sqrt(var(FP115L$TotalLogDistance1)/length(FP115L$TotalLogDistance1))
SEFP115L
mean(FP115L$TotalLogDistance1)

###SE for log Pesticide Mix 2 5ºC ###
FP25L<-filter(c, c$Trial== 'P Mix 2 5ºC')
SEFP25L<-sqrt(var(FP25L$TotalLogDistance1)/length(FP25L$TotalLogDistance1))
SEFP25L
mean(FP25L$TotalLogDistance1)

###SE for log Pesticide Mix 2 10ºC ###
FP210L<-filter(c, c$Trial== 'P Mix 2 10ºC ')
SEFP210L<-sqrt(var(FP210L$TotalLogDistance1)/length(FP210L$TotalLogDistance1))
SEFP210L
mean(FP210L$TotalLogDistance1)

###SE for log esticide Mix 2 15ºC ###
FP215L<-filter(c, c$Trial== 'P Mix 2 15ºC')
SEFP215L<-sqrt(var(FP215L$TotalLogDistance1)/length(FP215L$TotalLogDistance1))
SEFP215L
mean(FP215L$TotalLogDistance1)
#Statistical test??
##Linear Model###

table(d$Variance)


tempm<- lm(Distance~Temperature,data=d)
summary(tempm)

chemm<- lm(Distance~Chemical,data=d)
summary(chemm)

##Barchart for all (with chems in isolation) ##

c$cont<-c(rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Exp1", 12), rep("Exp1", 12),
          rep("Exp1", 12), rep("Exp2", 12),
          rep("Exp2", 12), rep("Exp2", 12),
          rep("Exp3", 12), rep("Exp3", 12),
          rep("Exp3", 12), rep("Exp4", 12),
          rep("Exp4", 12), rep("Exp4", 12),
          rep("Exp5", 12), rep("Exp5", 12),
          rep("Exp5", 12))
c$cont<-factor(c$cont, levels=c("Control", "Exp", "Exp1", "Exp2", "Exp3", "Exp4","Exp5"))
ggplot(c, aes(x=Trial, y=Distance.cm, fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC", "Fipronil 5ºC ", "Fipronil 10ºC", "Fipronil 15ºC ", "Fluralaner 5ºC ", "Fluralaner 10ºC ", "Fluralaner 15ºC ", "Selamectin 5ºC", "Selamectin 10ºC", "Selamectin 15ºC "
                                ,"P Mix 2 5ºC", "P Mix 2 10ºC ", "P Mix 2 15ºC","P Mix 1 5ºC ", "P Mix 1 10ºC", "P Mix 1 15ºC"))+
    scale_fill_manual(values=c("mediumorchid2", "turquoise","goldenrod2", "seagreen3", "blue", "purple", "red"))+
    theme(legend.position = "right", panel.background = element_blank())


##Just chems in isolation) ##

c$cont<-c(rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Exp1", 12), rep("Exp1", 12),
          rep("Exp1", 12), rep("Exp2", 12),
          rep("Exp2", 12), rep("Exp2", 12),
          rep("Exp3", 12), rep("Exp3", 12),
          rep("Exp3", 12), rep("Exp4", 12),
          rep("Exp4", 12), rep("Exp4", 12),
          rep("Exp5", 12), rep("Exp5", 12),
          rep("Exp5", 12))
c$cont<-factor(c$cont, levels=c("Control", "Exp", "Exp1", "Exp2", "Exp3", "Exp4","Exp5"))
ggplot(c, aes(x=Trial, y=Distance.cm, fill=cont)) +
    stat_summary(fun=mean, geom="bar")+
    stat_summary(fun.data=get.se, geom="errorbar", width=0.1)+
    scale_x_discrete(limits = c("Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC", "Fipronil 5ºC ", "Fipronil 10ºC", "Fipronil 15ºC ", "Fluralaner 5ºC ", 
                                "Fluralaner 10ºC ", "Fluralaner 15ºC ", "Selamectin 5ºC", "Selamectin 10ºC", "Selamectin 15ºC "))+
    scale_fill_manual(values=c("turquoise", "blue", "purple", "red"))+
    theme(legend.position = "none", panel.background = element_blank())





##Grouped by temp?##
c <- read.csv("Mobilitymaster.csv")

c$cont<-c(rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Exp1", 12), rep("Exp1", 12),
          rep("Exp1", 12), rep("Exp2", 12),
          rep("Exp2", 12), rep("Exp2", 12),
          rep("Exp3", 12), rep("Exp3", 12),
          rep("Exp3", 12), rep("Exp4", 12),
          rep("Exp4", 12), rep("Exp4", 12),
          rep("Exp5", 12), rep("Exp5", 12),
          rep("Exp5", 12))
c$cont<-factor(c$cont, levels=c("Control", "Exp", "Exp1", "Exp2", "Exp3", "Exp4","Exp5"))
ggplot(c, aes(x= Temperature, y=Distance.cm, fill=Chemical)) +
  geom_bar(position = "dodge", stat = "identity")+
    scale_x_discrete(limits = c("5ºC","10ºC","15ºC"))+ scale_fill_jama()
    

############### For errorbars ############
c <- read.csv("Mobilitymaster.csv")

# Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                       c(N    = length2(xx[[col]], na.rm=na.rm),
                         mean = mean   (xx[[col]], na.rm=na.rm),
                         sd   = sd     (xx[[col]], na.rm=na.rm)
                       )
                   },
                   measurevar
    )
    
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
}

library(plyr)
library(dplyr)

##### Does it work? ###
Ebars <- summarySE(c,measurevar= "Distance.cm", groupvars=c("Temperature", "Chemical"))
Ebars

##Grouped by chem?....... Works for error bars##  


Ebars2 <- Ebars
Ebars2$Chemical <- factor(Ebars2$Chemical)
ggplot(Ebars2, aes(x= Chemical, y=Distance.cm, fill = factor(Temperature, levels=c("5ºC", "10ºC", "15ºC")))) +
    geom_bar(position = "dodge", stat = "identity")+
    scale_x_discrete(limits = c("Control","Imidacloprid", "Selamectin", "Fluralaner", "Fipronil", "P Mix 2", "P Mix 1"))+
    labs(fill = "Temperature") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_errorbar(aes(ymin=Distance.cm-se, ymax=Distance.cm+se),  width=.22, position=position_dodge(.9)) +
  scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1')) +xlab("Temperature (ºC)") +
  ylab("Distance (cm)")



###Imidacloprid Concentrations works####

I <- read.csv("Imidacloprid Concentrations.csv")

ISUM <- summarySE(I,measurevar= "Distance.cm", groupvars=c("Temperature.ºC", "Concentration.µg.L"))
ISUM

ISUM2 <- ISUM 
ISUM2$Concentration.µg.L <- factor(ISUM2$Concentration.µg.L)
ggplot(ISUM2, aes(x= Concentration.µg.L, y=Distance.cm, fill = factor(Temperature.ºC, levels=c("5", "10", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("0","0.1","10","100"))+
  scale_fill_manual(values=c("steelblue", "lightskyblue1", "deepskyblue4")) +
  labs(fill = "Temperature.ºC", title = "Imidacloprid") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin=Distance.cm-se, ymax=Distance.cm+se),  width=.22, position=position_dodge(.9)) 


###P mIx all Concentrations works####

P <- read.csv("P Mix All Concentrations.csv")

PSUM <- summarySE(P,measurevar= "Distance.cm", groupvars=c("Temperature.ºC", "Concentration.µg.L"))
PSUM

PSUM2 <- PSUM 
PSUM2$Concentration.µg.L <- factor(PSUM2$Concentration.µg.L)
ggplot(PSUM2, aes(x= Concentration.µg.L, y=Distance.cm, fill = factor(Temperature.ºC, levels=c("5", "10", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("0","0.1","10","100"))+
  scale_fill_manual(values=c("steelblue", "lightskyblue1", "slateblue")) +
  labs(fill = "Temperature.ºC", title = "P Mix (all)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin=Distance.cm-se, ymax=Distance.cm+se),  width=.22, position=position_dodge(.9)) 




### Fipronil Concentrations Graph #####

Fip <- read.csv("Fipronil concentrations.xlsx.csv")

FIPSUM <- summarySE(Fip,measurevar= "Distance.cm", groupvars=c("Temperature.ºC", "Concentration.µg.L"))
FIPSUM

FIPSUM2 <- FIPSUM 
FIPSUM2$Concentration.µg.L <- factor(FIPSUM2$Concentration.µg.L)
ggplot(FIPSUM2, aes(x= Concentration.µg.L, y=Distance.cm, fill = factor(Temperature.ºC, levels=c("5", "10", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("0","0.1","10","100"))+
  scale_fill_manual(values=c("steelblue", "lightskyblue1", "slateblue")) +
  labs(fill = "Temperature.ºC", title = "Fipronil") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin=Distance.cm-se, ymax=Distance.cm+se),  width=.22, position=position_dodge(.9)) 

##Grouped by chem?##
c$cont<-c(rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Control", 12), rep("Exp", 12),
          rep("Exp1", 12), rep("Exp1", 12),
          rep("Exp1", 12), rep("Exp2", 12),
          rep("Exp2", 12), rep("Exp2", 12),
          rep("Exp3", 12), rep("Exp3", 12),
          rep("Exp3", 12), rep("Exp4", 12),
          rep("Exp4", 12), rep("Exp4", 12),
          rep("Exp5", 12), rep("Exp5", 12),
          rep("Exp5", 12))
c$cont<-factor(c$cont, levels=c("Control", "Exp", "Exp1", "Exp2", "Exp3", "Exp4","Exp5"))
ggplot(c, aes(x= Chemical, y=Distance.cm, fill = factor(Temperature, levels=c("5ºC", "10ºC", "15ºC")))) +
    geom_bar(position = "dodge", stat = "identity")+
    scale_x_discrete(limits = c("Control","Imidacloprid", "Selamectin", "Fluralaner", "Fipronil", "P Mix 2", "P Mix 1"))+
    scale_fill_manual(values=c("blue4", "lightskyblue1", "red3")) +
    labs(fill = "Temperature") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
                                                                                                                            



###Anova for all distance ##

aov <- aov(Distance.cm ~ Chemical * Temperature, data = c)
aov<- aov(Distance.cm ~ Chemical + Temperature + Chemical:Temperature, data = c)
summary(aov)

TukeyHSD(aov,which = "Chemical")
TukeyHSD(aov,conf.level=.95)

##Anova for all log distance ###

aov <- aov(TotalLogDistance1 ~ Chemical * Temperature, data = c)
aov<- aov(TotalLogDistance1 ~ Chemical + Temperature + Chemical:Temperature, data = c)
summary(aov)

TukeyHSD(aov,which = "Chemical")
TukeyHSD(aov,conf.level=.95)


