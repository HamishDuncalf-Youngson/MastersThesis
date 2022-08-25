
rm(list=ls())
install.packages('SciViews')
install.packages('wesanderson')
install.packages("gridExtra")
install.packages("tidyverse")
install.packages("ggsci")
install.packages("pals")

library(dplyr) 
library(SciViews)
library(ggplot2)
library(dplyr)
library(viridis)
library(wesanderson)
library(gridExtra)
library(tidyverse)
library(ggsci)
library(pals)



#Read in data
F <- read.csv('field bioassays2.csv')
Te <- read.csv('teabags.csv')
Ea <- read.csv('EA Concentraitions.csv')

#Coarse Leaf Litter k calculation 
F$coarseleafk1 <- (log(F$Post.Leaf.Litter.Coarse/F$Pre.Leaf.Litter.Weight..Coarse..g.)/21)*-1

#Fine Leaf Litter k calculation 
F$fineleafk1 <- (log(F$Post.Leaf.Litter.Fine/F$Pre.Leaf.Litter.Weight..Fine..g.)/21)*-1

#teabag 1 k calculation 
Te$teabagk <- (log(Te$Post.Teabag.Weight/Te$Pre.Teabag.Weight..g.)/21)*-1

#wettex fine k calculation (note the log +1!!!)
F$wettexfk <- (log(F$Post.Wettex..Fine/F$Pre.Wettex.Weight)/21)*-1

#wettex coarse k calculation (note the log +1!!!)
F$wettexck <- (log(F$Post.Wettex..Coarse/F$Pre.Wettex.Weight)/21)*-1

#######################Degree Days#############################################

#Degree Day Coarse Leaf Litter k calculation 
F$coarseleafddk1 <- (log(F$Post.Leaf.Litter.Coarse/F$Pre.Leaf.Litter.Weight..Coarse..g.)/(F$Temperature*21))*-1

#Degree Day Fine Leaf Litter k calculation 
F$fineleafddk1 <- (log(F$Post.Leaf.Litter.Fine/F$Pre.Leaf.Litter.Weight..Fine..g.)/(F$Temperature*21))*-1

#Degree Day teabag 1 k calculation 
Te$teabagddk <- (log(Te$Post.Teabag.Weight/Te$Pre.Teabag.Weight..g.)/(Te$Temperature*21))*-1

#Degree Day wettex fine k calculation (note the log +1!!!)
F$wettexfddk <- (log(F$Post.Wettex..Fine/F$Pre.Wettex.Weight)/(F$Temperature*21))*-1

#Degree Day wettex coarse k calculation (note the log +1!!!)
F$wettexcddk <- (log(F$Post.Wettex..Coarse/F$Pre.Wettex.Weight)/(F$Temperature*21))*-1

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
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



#Barchart for coarse leaf k 

leafcoarsesum<- summarySE(F,measurevar= "coarseleafk1", groupvars=c("Location"),na.rm = TRUE)

leafcoarsesum2 <- leafcoarsesum
  leafcoarsesum2$Location <- factor(leafcoarsesum2$Location)
coarseleaf <-ggplot(leafcoarsesum2, aes(x= Location, y=coarseleafk1, fill=Location)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Leaf Litter Decomposition Coarse") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= coarseleafk1-se, ymax= coarseleafk1+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Decay constant (k)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))

coarseleaf
#Barchart for fine leaf k 
leaffinesum<- summarySE(F,measurevar= "fineleafk1", groupvars=c("Location"),na.rm = TRUE)

leaffinesum2 <- leaffinesum
leaffinesum2$Location <- factor(leaffinesum2$Location)
fineleaf<- ggplot(leaffinesum2, aes(x= Location, y=fineleafk1, fill=Location)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Leaf Litter Decomposition Fine") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= fineleafk1-se, ymax= fineleafk1+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Decay constant (k)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))

fineleaf
#Barchart for Leaves  ####

leafcoarsesum<- summarySE(F,measurevar= "coarseleafk1", groupvars=c("Location"),na.rm = TRUE)
leafcoarsesum2 <- leafcoarsesum

leaffinesum<- summarySE(F,measurevar= "fineleafk1", groupvars=c("Location"),na.rm = TRUE)
leaffinesum2 <- leaffinesum

names(leafcoarsesum2)[3] <- "leafk" 
names(leaffinesum2)[3] <- "leafk" 
leafsum <-as.data.frame(rbind(leafcoarsesum2,leaffinesum2)) 
leafsum$Mesh <- c(rep("coarse", 3),rep("fine",3))

leafsum$Location <- factor(leafsum$Location)
Leaf <- ggplot(leafsum, aes(x= Location, y= leafk, fill=Location))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Leaf Litter Decomposition") + 
  theme(panel.background = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
        axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= leafk-se, ymax= leafk+se),  width=.22, position=position_dodge(.9)) +
  labs (x= "Location",  y="Decay constant (k)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1')) +
  facet_wrap(~Mesh)

Leaf



#Barchart for teabag k 
teabagsum1<- summarySE(Te,measurevar= "teabagk", groupvars=c("Location"),na.rm = TRUE)

teabagsum2 <- teabagsum1
teabagsum2$Location <- factor(teabagsum2$Location)
teabag<- ggplot(teabagsum2, aes(x= Location, y=teabagk, fill=Location)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Teabag Decomposition") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= teabagk-se, ymax= teabagk+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Decay constant (k)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))


teabag
#Barchart for fine wettex k 


wettexfinesum<- summarySE(F,measurevar= "wettexfk", groupvars=c("Location"),na.rm = TRUE)

wettexfinesum2 <- wettexfinesum
wettexfinesum2$Location <- factor(wettexfinesum2$Location)
finewettex <- ggplot(wettexfinesum2, aes(x= Location, y=wettexfk, fill=Location)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Fine Wettex Decomposition") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= wettexfk-se, ymax= wettexfk+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Decay constant (k)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))

finewettex

#Barchart for coarse wettex k
wettexcoarsesum<- summarySE(F,measurevar= "wettexck", groupvars=c("Location"),na.rm = TRUE)

wettexcoarsesum2 <- wettexcoarsesum
wettexcoarsesum2$Location <- factor(wettexcoarsesum2$Location)
coarsewettex<- ggplot(wettexcoarsesum2, aes(x= Location, y=wettexck, fill=Location)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Coarse Wettex Decomposition") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= wettexck-se, ymax= wettexck+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Decay constant (k)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))



#Barchart for Wettex  ####
wettexcoarsesum<- summarySE(F,measurevar= "wettexck", groupvars=c("Location"),na.rm = TRUE)
wettexcoarsesum2 <- wettexcoarsesum


wettexfinesum<- summarySE(F,measurevar= "wettexfk", groupvars=c("Location"),na.rm = TRUE)

wettexfinesum2 <- wettexfinesum


names(wettexfinesum2)[3] <- "wettexk" 
names(wettexcoarsesum2)[3] <- "wettexk" 
wettexsum <-as.data.frame(rbind(wettexcoarsesum2,wettexfinesum2)) 
wettexsum$Mesh <- c(rep("coarse", 3),rep("fine",3))

wettexfinesum <- wettexfinesum
wettexsum$Location <- factor(wettexsum$Location)
wettex <- ggplot(wettexsum, aes(x= Location, y= wettexk, fill=Location))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Wettex Decomposition") + 
  theme(panel.background = element_rect(fill = NA), 
    panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
                                    axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= wettexk-se, ymax= wettexk+se),  width=.22, position=position_dodge(.9)) +
  labs (x= "Location",  y="Decay constant (k)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1')) +
  facet_wrap(~Mesh)

wettex



##################################### DegreeDays ###############################

#Barchart for degree day coarse leaf k 
leafcoarseddsum<- summarySE(F,measurevar= "coarseleafddk1", groupvars=c("Location"),na.rm = TRUE)

leafcoarseddsum2 <- leafcoarseddsum
leafcoarseddsum2$Location <- factor(leafcoarseddsum2$Location)
degreedayscoarseleaf <- ggplot(leafcoarseddsum2, aes(x= Location, y=coarseleafddk1, fill=Location)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Leaf Litter Decomposition Coarse (Temperature Corrected)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= coarseleafddk1-se, ymax= coarseleafddk1+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Decay constant (kdd)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))

degreedayscoarseleaf

#Barchart for degree day fine leaf k 
leaffineddsum<- summarySE(F,measurevar= "fineleafddk1", groupvars=c("Location"),na.rm = TRUE)

leaffineddsum2 <- leaffineddsum
leaffineddsum2$Location <- factor(leaffineddsum2$Location)
degreedaysfineleaf <- ggplot(leaffineddsum2, aes(x= Location, y=fineleafddk1, fill=Location)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Leaf Litter Decomposition Fine (Temperature Corrected)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= fineleafddk1-se, ymax= fineleafddk1+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Decay constant (kdd)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))


degreedaysfineleaf

##Barchart for degree day leaves
leafcoarseddsum<- summarySE(F,measurevar= "coarseleafddk1", groupvars=c("Location"),na.rm = TRUE)
leafcoarseddsum2 <- leafcoarseddsum

leaffineddsum<- summarySE(F,measurevar= "fineleafddk1", groupvars=c("Location"),na.rm = TRUE)
leaffineddsum2 <- leaffineddsum


names(leafcoarseddsum2)[3] <- "leafk" 
names(leaffineddsum2)[3] <- "leafk" 
leafddsum <-as.data.frame(rbind(leafcoarseddsum2,leaffineddsum2)) 
leafddsum$Mesh <- c(rep("coarse", 3),rep("fine",3))

leafddsum$Location <- factor(leafddsum$Location)
Leafdd <- ggplot(leafddsum, aes(x= Location, y= leafk, fill=Location))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Leaf Litter Decomposition") + 
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
        axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= leafk-se, ymax= leafk+se),  width=.22, position=position_dodge(.9)) +
  labs (x= "Location",  y="Decay Constant (kDD-1)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1')) +
  facet_wrap(~Mesh)

Leafdd





#Barchart for degree day teabag k 
teabagddsum1<- summarySE(Te,measurevar= "teabagddk", groupvars=c("Location"),na.rm = TRUE)

teabagddsum2 <- teabagddsum1
teabagddsum2$Location <- factor(teabagddsum2$Location)
degreedaysteabag<- ggplot(teabagddsum2, aes(x= Location, y=teabagddk, fill=Location)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Teabag Decomposition ") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= teabagddk-se, ymax= teabagddk+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Decay constant (kDD-1)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))

degreedaysteabag


#Barchart for degree day fine wettex k 
wettexfineddsum<- summarySE(F,measurevar= "wettexfddk", groupvars=c("Location"),na.rm = TRUE)

wettexfineddsum2 <- wettexfineddsum
wettexfineddsum2$Location <- factor(wettexfineddsum2$Location)
degreedaysfinewettex<- ggplot(wettexfineddsum2, aes(x= Location, y=wettexfddk, fill=Location)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Fine Wettex Decomposition (Temperature Corrected)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= wettexfddk-se, ymax= wettexfddk+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Decay constant (kdd)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))

degreedaysfinewettex

#Barchart for degree day coarse wettex k 
wettexcoarseddsum<- summarySE(F,measurevar= "wettexcddk", groupvars=c("Location"),na.rm = TRUE)

wettexcoarseddsum2 <- wettexcoarseddsum
wettexcoarseddsum2$Location <- factor(wettexcoarseddsum2$Location)
degreedaycoursewettex <- ggplot(wettexcoarseddsum2, aes(x= Location, y=wettexcddk, fill=Location)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Coarse Wettex Decomposition (Temperature Corrected)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= wettexcddk-se, ymax= wettexcddk+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Decay constant (kdd)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))

degreedaycoursewettex


#Barchart for Wettex Degree days  ####
wettexcoarseddsum<- summarySE(F,measurevar= "wettexcddk", groupvars=c("Location"),na.rm = TRUE)
wettexcoarseddsum2 <- wettexcoarseddsum

wettexfineddsum<- summarySE(F,measurevar= "wettexfddk", groupvars=c("Location"),na.rm = TRUE)
wettexfineddsum2 <- wettexfineddsum

names(wettexfineddsum2)[3] <- "wettexk" 
names(wettexcoarseddsum2)[3] <- "wettexk" 
wettexddsum <-as.data.frame(rbind(wettexcoarseddsum2,wettexfineddsum2)) 
wettexddsum$Mesh <- c(rep("coarse", 3),rep("fine",3))


wettexddsum$Location <- factor(wettexddsum$Location)
wettexdd <- ggplot(wettexddsum, aes(x= Location, y= wettexk, fill=Location))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Wettex Decomposition") + 
  theme(panel.background = element_rect(fill = NA),  
        panel.grid.major = element_line(colour = "white"), axis.line = element_line(size = 0.5, colour = "grey"),
        axis.ticks = element_line(size = 0.5)) +
  geom_errorbar(aes(ymin= wettexk-se, ymax= wettexk+se),  width=.22, position=position_dodge(.9)) +
  labs (x= "Location",  y="Decay constant (kDD-1)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1')) +
  facet_wrap(~Mesh)

wettexdd





###degree days side- by-side with non #######
grid.arrange(wettex,wettexdd, ncol=2, nrow = 1)
grid.arrange(Leaf,Leafdd, ncol=2, nrow = 1)
grid.arrange(teabag,degreedaysteabag, ncol=2, nrow=1)


#Linear model for course leaf and temperature 

lc <- lm(F$coarseleafk1~F$Temperature)
summary(lc)
ggplot(F, aes(x=Temperature, y=coarseleafk1)) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Coarse Leaf Litter Decomposition", x= "Temperature (ºC)", y= "Decomposition (k)" ) 

lc <- lm(F$coarseleafk1~F$Temperature)
summary(lc)

#Scatterplot for Temperature and Fine Leaf Decomposition 

lf <- lm(F$fineleafk1~F$Temperature)
summary(lf)
ggplot(F, aes(x=Temperature, y=fineleafk1)) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Fine Leaf Litter Decomposition", x= "Temperature (ºC)", y= "Decomposition (k)" ) 

#Scatterplot for Temperature and Fine teabag decomposition

tb <- lm(Te$teabagk~Te$Temperature)
summary(tb)
ggplot(Te, aes(x=Temperature, y= teabagk)) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Teabag Decomposition", x= "Temperature (ºC)", y= "Decomposition (k)" ) 




#Scatterplot for Temperature and Fine Wettex Decomposition 

fw <- lm(F$wettexfk~F$Temperature)
summary(fw)
ggplot(F, aes(x=Temperature, y= wettexfk )) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Wettex Decomposition Fine", x= "Temperature (ºC)", y= "Decomposition (k)" ) 

#Scatterplot for Temperature and Coarse Wettex Decomposition 

cw <- lm(F$wettexck~F$Temperature)
summary(cw)
ggplot(F, aes(x=Temperature, y= wettexck )) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Wettex Decomposition Coarse", x= "Temperature (ºC)", y= "Decomposition (k)" ) 




####### By Location #########

#Boxplot for temperature by location
T <- read.csv('Average Temperature.csv')

ggplot(data=T, mapping=aes(x=Location, y=Temperature, fill= Location))+geom_boxplot() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values=c('turquoise','deepskyblue4', 'royalblue1')) + labs(title = "Temperature by location")


#Barchart for temperature by location
templocationsum<- summarySE(F,measurevar= "Temperature", groupvars=c("Location"),na.rm = TRUE)

templocationsum2 <- templocationsum 
templocationsum2$Location <- factor(templocationsum2$Location)
ggplot(templocationsum2, aes(x= Location, y=Temperature, fill=Location)) +   geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Temperature by location") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= Temperature-se, ymax= Temperature+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Temperature")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))

##### EA Concentration by location

##Boxplot for EA Concentration by location

ggplot(data=Ea, mapping=aes(x=Location, y=Imidacloprid.Concentration, fill= Location))+geom_boxplot() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values=c('turquoise','deepskyblue4', 'royalblue1')) + labs(title = "Imidacloprid")


#Barchart for EA imidacloprid concentration by location
imidalocationsum<- summarySE(Ea,measurevar= "Imidacloprid.Concentration", groupvars=c("Location"),na.rm = TRUE)

imidalocationsum2 <- imidalocationsum 
imidalocationsum2$Location <- factor(imidalocationsum2$Location)
Imidaclopridlocation <- ggplot(imidalocationsum2, aes(x= Location, y=Imidacloprid.Concentration, fill=Location)) +   geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Imidacloprid by location") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= Imidacloprid.Concentration-se, ymax= Imidacloprid.Concentration+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Imidacloprid Concentration (µg/L)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))

#Barchart for EA Fipronil concentration by location
fiplocationsum<- summarySE(Ea,measurevar= "Fipronil.Concentration", groupvars=c("Location"),na.rm = TRUE)

fiplocationsum2 <- fiplocationsum 
fiplocationsum2$Location <- factor(fiplocationsum2$Location)
Fipronillocation<- ggplot(fiplocationsum2, aes(x= Location, y=Fipronil.Concentration, fill=Location)) +   geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Fipronil by location") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= Fipronil.Concentration-se, ymax= Fipronil.Concentration+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Location",  y="Fipronil Concentration (µg/L)")  + scale_fill_manual(values = c('turquoise','deepskyblue4', 'royalblue1'))

Fipronillocation



######################

#Scatterplot for Gammarus and Coarse Leaf Decomposition 

gcl <- lm(F$coarseleafddk1~F$Gammarus)
summary(gcl)
ggplot(F, aes(x=Gammarus, y=coarseleafddk1 )) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Leaf Decomposition Coarse Gammarus", x= "Number of Gammarus", y= "Decomposition (kdd)" ) 



#Scatterplot for Gammarus and Coarse Wettex Decomposition
gcw <- lm(F$wettexfddk~F$Gammarus)
summary(gcw)
ggplot(F, aes(x=Gammarus, y=wettexcddk )) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Wettex Decomposition Coarse Gammarus", x= "Number of Gammarus", y= "Decomposition (kdd)" ) 



######### Using EA chemical data ##############################################

####Coarse Wettex

##Imidacloprid Coarse Wettex 
icw <- lm(F$wettexcddk~F$Imidacloprid.Concentration)
summary(icw)
ggplot(F, aes(x=Imidacloprid.Concentration, y=wettexcddk )) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Wettex Decomposition Coarse Imidacloprid", x= "Imidacloprid Concentration (µg/L)", y= "Decomposition (kdd)" ) 

imidaclopridcwettex
## Fipronil Coarse Wettex 

fcw <- lm(F$wettexcddk~F$Fipronil.Concentration)
summary(fcw)
ggplot(F, aes(x=Fipronil.Concentration, y=wettexcddk )) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Wettex Decomposition Coarse Fipronil", x= "Fipronil Concentration (µg/L)", y= "Decomposition (kdd)" ) 

fipronilcwettex

#### Fine Wettex


#Imidacloprid Fine Wettex 
ifw <- lm(F$wettexfddk~F$Imidacloprid.Concentration)
summary(ifw)
 ggplot(F, aes(x=Imidacloprid.Concentration, y=wettexfddk )) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Wettex Decomposition Fine Imidacloprid", x= "Imidacloprid Concentration (µg/L)", y= "Decomposition (kdd)" ) 


#Fipronil Fine Wettex 
ffw <- lm(F$wettexfddk~F$Fipronil.Concentration)
summary(ffw)
ggplot(F, aes(x=Fipronil.Concentration, y=wettexfddk )) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Wettex Decomposition Fine Fipronil", x= "Fipronil Concentration (µg/L)", y= "Decomposition (kdd)" ) 


#### Coarse Leaf

#Imidacloprid Coarse Leaf 
icl <- lm(F$coarseleafddk1~F$Imidacloprid.Concentration)
summary(icl)
ggplot(F, aes(x=Imidacloprid.Concentration, y=coarseleafddk1)) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Leaf Litter Decomposition Coarse Imidacloprid", x= "Imidacloprid Concentration (µg/L)", y= "Decomposition (kdd)" ) 

#Fipronil Coarse Leaf 
fcl <- lm(F$coarseleafddk1~F$Fipronil.Concentration)
summary(fcl)
ggplot(F, aes(x=Fipronil.Concentration, y=coarseleafk1 )) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Leaf Litter Decomposition Coarse Fipronil", x= "Fipronil Concentration (µg/L)", y= "Decomposition (kdd)" ) 


### Fine Leaf

#Imidacloprid Fine Leaf
ifl <- lm(F$fineleafddk1~F$Imidacloprid.Concentration)
summary(ifl)
ggplot(F, aes(x=Imidacloprid.Concentration, y=fineleafddk1 )) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Leaf Litter Decomposition Fine Imidacloprid", x= "Imidacloprid Concentration (µg/L)", y= "Decomposition (kdd)" ) 


#Fipronil Fine Leaf 
ffl <- lm(F$fineleafddk1~F$Fipronil.Concentration)
summary(ffl)
ggplot(F, aes(x=Imidacloprid.Concentration, y=fineleafddk1 )) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Leaf Litter Decomposition Fine Fipronil", x= "Fipronil Concentration (µg/L)", y= "Decomposition (kdd)" ) 


##### Teabag

#Imidacloprid Teabag 
ite <- lm(Te$teabagddk~Te$Imidacloprid.Concentration)
summary(ite)
ggplot(Te, aes(x=Imidacloprid.Concentration, y=teabagddk )) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Teabag Decomposition Imidacloprid", x= "Imidacloprid Concentration (µg/L)", y= "Decomposition (kdd)" ) 

#Fipronil Teabag 
fte <- lm(Te$teabagddk~Te$Fipronil.Concentration)
summary(fte)
ggplot(Te, aes(x=Fipronil.Concentration, y=teabagddk )) + geom_point(size=2.5, shape=18) + geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Teabag Decomposition Fipronil", x= "Fipronil Concentration (µg/L)", y= "Decomposition (kdd)" ) 

######


#Barchart for degree day coarse leaf k  by site?? 
leafcoarseddsumriv<- summarySE(F,measurevar= "coarseleafddk1", groupvars=c("Site"),na.rm = TRUE)
leafcoarseddsumriv$Location <- c(rep("Norfolk",2),rep("London",3), rep("Norfolk"), rep("London",2), 
rep("Wiltshire",4), rep("Norfolk", 3), rep("Wiltshire",2), rep("London"))


leafcoarseddsumriv2 <- leafcoarseddsumriv
leafcoarseddsumriv2$Site <- factor(leafcoarseddsumriv2$Site, levels = c("Crane","Cray","Darent","Hogsmill","Lea", "Wandle", 
                                                                        "Bayfield Hall", "Bintry Mill Trout Fishery", "Glandford", "Pensthorpe", "Sculthorpe","Selbrigg", "Littlecote", "Manton Bridge", "Marlborough College", "Mildenhall Trout Farm", 
                                                                        "Stitchcombe", "Stonebridge Meadow"))

degreedayscoarseleafriv <- ggplot(leafcoarseddsumriv2, aes(x= Site, y=coarseleafddk1, fill=Site)) +
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "Leaf Litter Decomposition Coarse") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= coarseleafddk1-se, ymax= coarseleafddk1+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Site",  y="Decay constant (kDD-1)") + facet_wrap(~Location, scales="free")  + 
  scale_fill_manual(values = c(rep('turquoise',6), rep('deepskyblue4',6), rep('royalblue1',6))) +
  theme(axis.text.x = element_text(angle = 345))

degreedayscoarseleafriv






