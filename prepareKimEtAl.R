

library(reshape2)


m.kim <- read.csv(file = './externaldata/KimCell08/mmc2.csv')
m.kim <- melt(m.kim, c('Symbol','Accession.number')) 
m.kim <- subset(m.kim, value == 'target')
m.kim$meta <- m.kim$variable
m.kim$variable <- 'GeneSymbol'
m.kim$source <- 'kimTF'
m.kim$value <- m.kim$Symbol
m.kim <- m.kim[,c('source','meta','variable','value')]
save( list = 'm.kim', file = './intermediate_results/m.kim.RData')




