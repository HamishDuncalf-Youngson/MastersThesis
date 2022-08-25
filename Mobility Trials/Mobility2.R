rm(list=ls())


c <- read.csv("Mobilitymaster.csv")
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
distancetemperature<- summarySE(c,measurevar= "Distance.cm", groupvars=c("Temperature"))
distancetemperature2 <- distancetemperature
distancetemperature2$Temperature <- factor(distancetemperature2$Temperature,levels = c("5ºC", "10ºC", "15ºC"))
ggplot(distancetemperature2, aes(x= Temperature, y=Distance.cm, fill=Temperature)) +   geom_bar(position = "dodge", stat = "identity") + 
  geom_errorbar(aes(ymin= Distance.cm-se, ymax= Distance.cm+se),  width=.22, position=position_dodge(.9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs (x= "Temperature (ºC)",  y="Distance (cm)")  + scale_fill_manual(values = c('deepskyblue4','steelblue', 'lightskyblue1'))




#Barchart for all Distance by chemical
distancechemical<- summarySE(c,measurevar= "Distance.cm", groupvars=c("Chemical"),na.rm = TRUE)
distancechemical2 <- distancechemical
distancechemical2$Chemical <- factor(distancechemical2$Chemical)
ggplot(distancechemical2, aes(x= Chemical, y=Distance.cm, fill=Chemical)) +   geom_bar(position = "dodge", stat = "identity") + 
  geom_errorbar(aes(ymin= Distance.cm-se, ymax= Distance.cm+se),  width=.22, position=position_dodge(.9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs (x= "Chemical Treatment",  y="Distance (cm)")  + scale_fill_lancet()


###Barchart for total distance Imidacloprid###
###and whichever other ones I want to compare to control.... Just pop in whichever chemical you need to look at. 

c <- read.csv("Mobilitymaster.csv")

c$cont<-c(rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6),
          rep("Control", 6), rep("Exp", 6))
c$cont<-factor(c$cont, levels=c("Control", "Exp"))

trialdistance<- summarySE(c,measurevar= "Distance.cm", groupvars=c("Trial"), na.rm = TRUE)
trialdistance2 <- trialdistance
trialdistance2$Trial <- factor(trialdistance2$Trial)

trialdistance2$chemical <- c(rep("Control", 3),rep("Fipronil",3),rep("Fluralaner",3), rep("Imidacloprid",3), rep("P Mix 1", 3), rep("P Mix 2", 3), rep("Selamectin",3))

ggplot(trialdistance2, aes(x=Trial, y=Distance.cm, fill=chemical)) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(title = "Imidacloprid 100 (µg/L)", x= "",  y="Distance (cm)") +
  geom_errorbar(aes(ymin= Distance.cm-se, ymax= Distance.cm+se),  width=.22, position=position_dodge(.9)) +
  scale_x_discrete(limits = c("Control 5ºC","Control 10ºC", "Control 15ºC","Imidacloprid 5ºC","Imidacloprid 10ºC","Imidacloprid 15ºC"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none")



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



### All chems one graph

get.se <- function(y) {
  se <- sd(y)/sqrt(length(y))
  mu <- mean(y)
  c(ymin=mu-se, ymax=mu+se)
}

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
  theme(legend.position = "", panel.background = element_blank(), axis.text.x= element_text(angle=340))  + scale_fill_jama()


##Grouped by temp? 


triallyy2 <- summarySE(c,measurevar= "Distance.cm", groupvars=c("Temperature", "Chemical"),na.rm = TRUE)
triallyy2$Temperature <- factor(triallyy2$Temperature)

ggplot(triallyy2, aes(x= Temperature, y=Distance.cm, fill=Chemical)) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("5ºC","10ºC","15ºC"))+ scale_fill_jama() +theme(panel.background = element_blank()) +
  geom_errorbar(aes(ymin= Distance.cm-se, ymax= Distance.cm+se),width=0.22, position=position_dodge(.9)) +xlab("Temperature (ºC)") +
  ylab("Distance (cm)")





##Grouped by chem?

Ebars <- summarySE(c,measurevar= "Distance.cm", groupvars=c("Temperature", "Chemical"))
Ebars

Ebars2 <- Ebars
Ebars2$Chemical <- factor(Ebars2$Chemical)
ggplot(Ebars2, aes(x= Chemical, y=Distance.cm, fill = factor(Temperature, levels=c("5ºC", "10ºC", "15ºC")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("Control","Imidacloprid", "Selamectin", "Fluralaner", "Fipronil", "P Mix 2", "P Mix 1"))+
  labs(fill = "Temperature") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin=Distance.cm-se, ymax=Distance.cm+se),  width=.22, position=position_dodge(.9)) +
  scale_fill_manual(values = c('deepskyblue4','steelblue', 'lightskyblue1')) +xlab("Chemical") +
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
  scale_fill_manual(values=c('deepskyblue4','steelblue', 'lightskyblue1')) +
  labs(fill = "Temperature.ºC", title = "Imidacloprid") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin=Distance.cm-se, ymax=Distance.cm+se),  width=.22, position=position_dodge(.9)) +xlab("Concentration (µg/L)") +
  ylab("Distance (cm)")


###P mix all Concentrations works####

P <- read.csv("P Mix All Concentrations.csv")

PSUM <- summarySE(P,measurevar= "Distance.cm", groupvars=c("Temperature.ºC", "Concentration.µg.L"))
PSUM

PSUM2 <- PSUM 
PSUM2$Concentration.µg.L <- factor(PSUM2$Concentration.µg.L)
ggplot(PSUM2, aes(x= Concentration.µg.L, y=Distance.cm, fill = factor(Temperature.ºC, levels=c("5", "10", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("0","0.1","10","100"))+
  scale_fill_manual(values=c('deepskyblue4','steelblue', 'lightskyblue1')) +
  labs(fill = "Temperature.ºC", title = "P Mix (all)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin=Distance.cm-se, ymax=Distance.cm+se),  width=.22, position=position_dodge(.9)) +xlab("Concentration (µg/L)") +
  ylab("Distance (cm)")


### Fipronil Concentrations Graph #####

Fip <- read.csv("Fipronil concentrations.xlsx.csv")

FIPSUM <- summarySE(Fip,measurevar= "Distance.cm", groupvars=c("Temperature.ºC", "Concentration.µg.L"))
FIPSUM

FIPSUM2 <- FIPSUM 
FIPSUM2$Concentration.µg.L <- factor(FIPSUM2$Concentration.µg.L)
ggplot(FIPSUM2, aes(x= Concentration.µg.L, y=Distance.cm, fill = factor(Temperature.ºC, levels=c("5", "10", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("0","0.1","10","100"))+
  scale_fill_manual(values=c('deepskyblue4','steelblue', 'lightskyblue1')) +
  labs(fill = "Temperature.ºC", title = "Fipronil") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin=Distance.cm-se, ymax=Distance.cm+se),  width=.22, position=position_dodge(.9)) +xlab("Concentration (µg/L)") +
  ylab("Distance (cm)")







