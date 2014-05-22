dir.create(file.path('outputdata/figReconstitution'), showWarnings = FALSE)

library(ggplot2)
library(plyr)
library(reshape2)
library(stringr)



reconstitutionData <- read.csv('./inputdata/FACSreconstitution.csv')
dat <- ddply(reconstitutionData, c('Experiment', 'time'), summarize, 
    avgFracPLUS = mean(fracPLUS, na.rm = TRUE),
    stdFracPLUS = sd(fracPLUS, na.rm = TRUE))


# \section{Stabilization in 2i/LIF}


dftoplot <- subset(dat, Experiment == 'fromSerumLIF')
q <- ggplot(dftoplot, aes(x = time)) +
    geom_line(aes(y = avgFracPLUS, group = Experiment)) + 
    geom_point(aes(y = avgFracPLUS)) +
    geom_errorbar(aes(ymin = avgFracPLUS - stdFracPLUS, 
	ymax = avgFracPLUS + stdFracPLUS), width=3) +
    scale_x_continuous(breaks = seq(from=0, to=max(dat$time,na.rm=TRUE), by=12)) +
    ylim(0,100) + ylab('% VNP(+) cells')
q <- q + xlab('hours after transfer to 2i+LIF')
q <- q + theme_minimal()



setEPS()
postscript('./outputdata/figReconstitution/transferTo2iLIF.eps', width=2.8, height=2.4)
print(q)
dev.off()


# \section{Reconstitution from purified populations}


dftoplot <- subset(dat, Experiment %in% c('fromGFP','fromNEG'))
q <- ggplot(dftoplot, aes(x = time)) +
    geom_line(aes(y = avgFracPLUS, group = Experiment, color = Experiment)) +
    geom_errorbar(aes(ymin = avgFracPLUS - stdFracPLUS, 
	ymax = avgFracPLUS + stdFracPLUS), width=3) +
    geom_point(aes(y = avgFracPLUS)) +
    scale_x_continuous(breaks = seq(from=0, to=max(dat$time,na.rm=TRUE), by=24)) +
    ylim(0,100) + ylab('% VNP(+) cells') + scale_color_discrete(guide = FALSE)
q <- q + xlab('hours after sorting')
q <- q + theme_minimal()
setEPS()
postscript('./outputdata/figReconstitution/reconstitution.eps', width=2.8, height=2.4)
print(q)
dev.off()

