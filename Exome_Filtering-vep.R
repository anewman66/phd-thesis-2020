setwd("/users/alex/Thesis_Writing/Exome/PrognosticAnalysis")

M3q29 <- read.table(file="3q29_mutations_39_cases_paired_format.txt",sep='\t',header=T)
M17q <- read.table(file="17q_mutations_39_cases_paired_format.txt",sep='\t',header=T)

U3q29 <- read.table(file="3q29_mutations_39_cases_unpaired.txt",sep='\t',header=T)
U17q <- read.table(file="17q_mutations_39_cases_unpaired.txt",sep='\t',header=T)


table(M3q29$SYMBOL)
sort(table(M3q29$SYMBOL))

table(M17q$SYMBOL)
sort(table(M17q$SYMBOL))

sort(table(M17q$Sample))
unique(table(M17q$Sample))
length(unique(table(M17q$Sample)))
dim(M17q)

M17q.filt1 <- M17q[M17q$IMPACT != "MODIFIER",]
M17q.filt2 <- M17q.filt1[M17q.filt1$IMPACT  != "LOW",]

sort(table(M17q.filt2$Sample))
unique(table(M17q.filt2$Sample))
length(unique(table(M17q.filt2$Sample)))
dim(M17q.filt2)
sort(table(M17q.filt2$SYMBOL))
length(sort(table(M17q.filt2$SYMBOL)))
M17q.genes <- sort(table(M17q.filt2$SYMBOL))
M17q.genes1 <- M17q.genes[M17q.genes != 0]
M17q.genes2 <- M17q.genes1[M17q.genes1 != 1]

c <-M17q.filt2[M17q.filt2$SYMBOL == "RNF213",]
c$Sample

#Unpaired 17q
sort(table(U17q$SYMBOL))

sort(table(U17q$Sample))
unique(table(U17q$Sample))
length(unique(table(U17q$Sample)))
dim(U17q)

U17q.filt1 <- U17q[U17q$IMPACT != "MODIFIER",]
U17q.filt2 <- U17q.filt1[U17q.filt1$IMPACT  != "LOW",]

sort(table(U17q.filt2$Sample))
unique(table(U17q.filt2$Sample))
length(unique(table(U17q.filt2$Sample)))
dim(U17q.filt2)
sort(table(U17q.filt2$SYMBOL))
length(sort(table(U17q.filt2$SYMBOL)))
U17q.genes <- sort(table(U17q.filt2$SYMBOL))
U17q.genes1 <- U17q.genes[U17q.genes != 0]
U17q.genes2 <- U17q.genes1[U17q.genes1 != 1]

ca <-U17q.filt2[U17q.filt2$SYMBOL == "RNF213",]
ca$Sample

#3q29
sort(table(M3q29$Sample))
unique(table(M3q29$Sample))
length(unique(table(M3q29$Sample)))
dim(M3q29)

M3q29.filt1 <- M3q29[M3q29$IMPACT != "MODIFIER",]
M3q29.filt2 <- M3q29.filt1[M3q29.filt1$IMPACT  != "LOW",]

dim(M3q29.filt2)

sort(table(M3q29.filt2$Sample))
unique(table(M3q29.filt2$Sample))
length(unique(table(M3q29.filt2$Sample)))
dim(M3q29.filt2)
sort(table(M3q29.filt2$SYMBOL))
length(sort(table(M3q29.filt2$SYMBOL)))
M3q29.genes <- sort(table(M3q29.filt2$SYMBOL))
M3q29.genes1 <- M3q29.genes[M3q29.genes != 0]
M3q29.genes2 <- M3q29.genes1[M3q29.genes1 != 1]

c <-M3q29.filt2[M3q29.filt2$SYMBOL == "",]
c$Sample



#Unpaired 17q
sort(table(U3q29$SYMBOL))

sort(table(U3q29$Sample))
unique(table(U3q29$Sample))
length(unique(table(U3q29$Sample)))
dim(U3q29)

U3q29.filt1 <- U3q29[U3q29$IMPACT != "MODIFIER",]
U3q29.filt2 <- U3q29.filt1[U3q29.filt1$IMPACT  != "LOW",]

sort(table(U3q29.filt2$Sample))
unique(table(U3q29.filt2$Sample))
length(unique(table(U3q29.filt2$Sample)))
dim(U3q29.filt2)
sort(table(U3q29.filt2$SYMBOL))
length(sort(table(U3q29.filt2$SYMBOL)))
U3q29.genes <- sort(table(U3q29.filt2$SYMBOL))
U3q29.genes1 <- U3q29.genes[U3q29.genes != 0]
U3q29.genes2 <- U3q29.genes1[U3q29.genes1 != 1]

ca <-U3q29.filt2[U3q29.filt2$SYMBOL == "MUC4",]
U3q29.filt3 <- U3q29.filt2[U3q29.filt2$SYMBOL == "MUC4",]
length(sort(table(U3q29.filt3$Existing_variation)))
sort(table(U3q29.filt3$Existing_variation))

ca$Sample

