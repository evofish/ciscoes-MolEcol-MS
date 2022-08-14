library(dplyr)

setwd('~/Desktop/')
x<-read.table('go-test.txt', sep='\t', header=T) 

head(x) #make sure the table has headers 

go<-x %>% group_by(gene) %>% summarize(text= paste(GO, collapse=";"))

write.table(go, file="summarized-GO.txt", sep='\t')
