rm(list=ls())

library(plyr)
library(dplyr)
library(SciViews)
library(ggplot2)
library(dplyr)

#Calculating k for leaf litter weight with Gammarus
L <- read.csv("Leaf litter weights2 .csv")

kwlog<- log(L$Weight..g)
kwlog

LK <- mutate(L,kdecay=kwlog/21)
LK

LK1 <- mutate(LK,kdecay1=kdecay*-1)
LK1


#Calculating k for leaf litter weight without Gammarus
LC <- read.csv("Leaf Litter sans inverts.csv")

LCkwlog<- log(LC$Weight..g)
LCkwlog

LCLK <- mutate(LC,kdecay=LCkwlog/21)
LCLK

LCLK1 <- mutate(LCLK,kdecay1=kdecay*-1)
LCLK1

#Calculating k for teabag weight with Gammarus
T <- read.csv("Tea bag weights3.csv")

tkwlog<- log(T$Weight..g./1.9)
tkwlog

TLK <- mutate(T,kdecay=tkwlog/21)
TLK

TLK1 <- mutate(TLK,kdecay1=kdecay*-1)
TLK1

#Calculating k for teabag weight without Gammarus
TC <- read.csv("Tea bags sans inverts.csv")

TCtkwlog<- log(TC$Weight..g./1.9)
TCtkwlog

TCTLK <- mutate(TC,kdecay=TCtkwlog/21)
TCTLK

TCTLK1 <- mutate(TCTLK,kdecay1=kdecay*-1)
TCTLK1


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
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


#Barchart for k leaf litter with Gammarus

LK1sum<- summarySE(LK1,measurevar= "kdecay1", groupvars=c("Temperature..CºC.", "Treatment"))

LK1sum2 <- LK1sum
LK1sum2$Treatment <- factor(LK1sum2$Treatment)
ggplot(LK1sum2, aes(x= Treatment, y= kdecay1, fill = factor(Temperature..CºC., levels=c("5", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("Control", "Imidacloprid", "Fipronil", "P Mix"))+
  scale_fill_manual(values=c("deepskyblue4", "lightskyblue1", "slateblue","lightskyblue1")) +
  labs(fill = "Temperature ºC", title = "Leaf Litter Trials - Gammarus") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= kdecay1-se, ymax= kdecay1+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Chemical",  y= expression(paste("Decay constant k " (DD^-1)))) +
  ylim(0, 0.02)


#Barchart for k leaf litter without Gammarus

LCLK1sum<- summarySE(LCLK1,measurevar= "kdecay1", groupvars=c("Temperature..CºC.", "Treatment"))

LCLK1sum2 <- LCLK1sum
LCLK1sum2$Treatment <- factor(LCLK1sum2$Treatment)
ggplot(LCLK1sum2, aes(x= Treatment, y= kdecay1, fill = factor(Temperature..CºC., levels=c("5", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("Control", "Imidacloprid", "Fipronil", "P Mix"))+
  scale_fill_manual(values=c("deepskyblue4", "lightskyblue1")) +
  labs(fill = "Temperature ºC", title = "Leaf Litter Trials - No Gammarus") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= kdecay1-se, ymax= kdecay1+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Chemical", y= expression(paste("Decay constant k " (DD^-1))))+ ylim(0, 0.02)



#Barchart for k tea bags with Gammarus

TLK1sum<- summarySE(TLK1,measurevar= "kdecay1", groupvars=c("Temperature..ºC.", "Treatment"))

TLK1sum2 <- TLK1sum
TLK1sum2$Treatment <- factor(TLK1sum2$Treatment)
ggplot(TLK1sum2, aes(x= Treatment, y= kdecay1, fill = factor(Temperature..ºC., levels=c("5", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("Control", "Imidacloprid","Fipronil", "P Mix"))+
  scale_fill_manual(values=c("deepskyblue4", "lightskyblue1", "slateblue","skylightblue1")) + 
  labs(fill = "Temperature ºC", title = "Green Teabags - Gammarus") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= kdecay1-se, ymax= kdecay1+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Chemical", y= expression(paste("Decay constant k " (DD^-1))))+
  ylim(0, 0.04)



#Barchart for k tea bags without Gammarus

TCTLK1sum<- summarySE(TCTLK1,measurevar= "kdecay1", groupvars=c("Temperature..ºC.", "Treatment"))

TCTLK1sum2 <- TCTLK1sum
TCTLK1sum2$Treatment <- factor(TCTLK1sum2$Treatment)
ggplot(TCTLK1sum2, aes(x= Treatment, y= kdecay1, fill = factor(Temperature..ºC., levels=c("5", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("Control", "Imidacloprid", "Fipronil", "P Mix"))+
  scale_fill_manual(values=c("deepskyblue4", "lightskyblue1", "steelblue","lightskyblue")) + 
  labs(fill = "Temperature ºC", title = "Green Teabags - No Gammarus") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin= kdecay1-se, ymax= kdecay1+se),  width=.22, position=position_dodge(.9))  +
  labs (x= "Chemical", y= expression(paste("Decay constant k " (DD^-1))))+
  ylim(0, 0.04)



