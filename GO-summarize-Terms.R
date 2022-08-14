##script for summarizing the GO annotation for the C. artedi transcriptome. 
##the script for the GO analysis is available at: https://github.com/z0on/GO_MWU

library(dplyr)

setwd('~/Desktop/')
x<-read.table('go-test.txt', sep='\t', header=T) 

head(x) #make sure the table has headers 

go<-x %>% group_by(gene) %>% summarize(text= paste(GO, collapse=";"))

write.table(go, file="summarized-GO.txt", sep='\t')
