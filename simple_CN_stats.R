setwd('/Users/alex/Thesis_Writing/Chapter6_eBL/eBL_GISTIC/')
dex <- read.table('All_sample_LOH_seg.txt',sep='\t',header=T)

str(dex)
dex$size <- dex$End.Position - dex$Start.Position
dex2 <- dex[dex$size > 1000000,]
dex3a <- dex2[dex2$Seg.CN > 0.15,]
dex3b <- dex2[dex2$Seg.CN < -0.15,]

dex3 <- rbind(dex3a,dex3b)

cn <-table(dex3$Sample)
median(cn)
mean(cn)
range(cn)
