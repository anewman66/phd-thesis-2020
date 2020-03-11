#Forest plot tests

BiocManager::install('forestplot')
library(forestplot)

setwd('/Users/alex/Thesis_Writing/Chapter7_Pairs')

input <- read.table('forest_input.txt',sep = '\t', header = T)

#With no table
forestplot(c(1,2,3,4,5,6,7,8,9,10,11),
           input$Relapse,
           input$Diagnosis,
           input$LFU,
           clip=c(0,36),
           xlab="Time (months)",
           xticks=c(0,6,12,18,24,30,36),
           boxsize = 0.1)

#With table
forestplot(input[,-3],
           mean=input$Relapse,
           lower=input$Diagnosis,
           upper=input$LFU,
           clip=c(0,36),
           xlab="Time (months)",
           xticks=c(0,6,12,18,24,30,36),
           boxsize = 0.1,
           txt_gp=fpTxtGp(label=gpar(fontfamily='',cex=1)),
           align='l')


?forestplot
