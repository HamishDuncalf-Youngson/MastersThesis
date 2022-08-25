
rm(list=ls())


P <-read.csv("PredatorPrey Trials 2.csv")

library(ggplot2)
library(dplyr)


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


P$Feedingrate <- P$Prey.eaten/24

Psum <- summarySE(P,measurevar= "Feedingrate", groupvars=c("Temperature", "Treatment"))
Psum

Psum2 <- Psum
Psum2$Treatment <- factor(Psum2$Treatment)
ggplot(Psum2, aes(x= Treatment, y= Feedingrate, fill = factor(Temperature, levels=c("5", "15")))) +
  geom_bar(position = "dodge", stat = "identity")+
  scale_x_discrete(limits = c("Control", "Imidacloprid", "Fipronil", "P Mix"))+
  scale_fill_manual(values=c("deepskyblue4", "lightskyblue1", "slateblue")) +
  labs(fill = "Temperature ºC", title = "Predator-Prey Trials") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin=Feedingrate-se, ymax=Feedingrate+se),  width=.22, position=position_dodge(.9)) + labs (x= "",  y="Feeding Rate (individuals per hour)") 






