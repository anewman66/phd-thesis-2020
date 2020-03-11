source("https://bioconductor.org/biocLite.R")

#You might not need to install any of these other than GenVisR, I had issues with dependencies so had to install some packages separately.

BiocManager::install(c('GenVisR'))

biocLite("GenVisR")

library('GenVisR')

#Set working directory
setwd('/Users/alex/Genomics/Malawi/Lolliplot')

#Read in the file with mutation information. Make sure each mutation is on a separate line.
sBL.muts <- read.table(file="sBL_muts.txt",sep='\t',header=T)
#I think it already creates a dataframe but it seemed to work more often if I did this.
sBL.muts <- as.data.frame(sBL.muts)

eBL.muts <- read.table(file="eBL_muts.txt",sep='\t',header=T)
#I think it already creates a dataframe but it seemed to work more often if I did this.
eBL.muts <- as.data.frame(eBL.muts)

#all.muts <- read.table(file='all_TP53_muts.txt',sep='\t',header=T)

#Provide a file of domains for your gene. You can get this from ensembl, but it will likely have more domains than you need, so edit it.
tp53.dom <- read.table(file="TP53_ensembl_domains.txt",sep='\t',header=T)
#Again, might just make it run a bit better.
tp53.dom <- as.data.frame(tp53.dom)

#For one group, coloured by subtype
#lolliplot(CCLG.muts, labelCol = "amino_acid_change",fillCol="Diagnosis",txtSize = 4, txtAngle = 90,z=tp53.dom)

#For MYC+ on top, MYC- and unk underneath
lolliplot(eBL.muts,sBL.muts, labelCol = "amino_acid_change",txtSize = 4, txtAngle = 90,z=tp53.dom)

#This makes the actual plot.
#X is the mutation table.
#Y is not used here, but can be used to plot more mutations underneath the plot.
#Z is the domain file.
#labelCol will label each mutation with a column from your input file. I chose amino acid change. 
#txtSize is just there to play around with, same as angle. I'd delete them for your first plot.

lolliplot(x= all.muts, labelCol = "amino_acid_change",txtSize = 4, txtAngle = 90,z=tp53.dom)
#pdf(file='demo_file.pdf',plot)
#The plot does take several minutes to plot, and will sometimes come back with a multitude of errors!
#Re-run it without changing anything and eventually it will work again. I don't know why but it works.
#Why does it run so much better on the mac?



??GenVisR
??lolliplot

#FOXO1

setwd('N:/Alex Newman/R_workspace/GenVisR')
foxo1.muts1 <- read.table(file="FOXO1_literature_plot.txt",sep='\t',header=T)
foxo1.muts1 <- as.data.frame(foxo1.muts1)

foxo1.muts.eBL <- read.table(file="FOXO1_endemic_BL_plot.txt",sep='\t',header=T)
foxo1.muts.eBL <- as.data.frame(foxo1.muts.eBL)

foxo1.muts.sBL <- read.table(file="FOXO1_sporadic_BL_plot.txt",sep='\t',header=T)
foxo1.muts.sBL <- as.data.frame(foxo1.muts.sBL)

foxo1.dom <- read.table(file="FOXO1_domains.txt",sep='\t',header=T)
foxo1.dom <- as.data.frame(foxo1.dom)
foxo1.dom2 <- read.table(file="FOXO1_domains_exons.txt.",sep='\t',header=T)

#Lit plot A
lolliplot(CCLG.muts, labelCol = "amino_acid_change",fillCol="Diagnosis",txtSize = 4, txtAngle = 90,z=tp53.dom)


#Our plot B
lolliplot(foxo1.muts.sBL,foxo1.muts.eBL, labelCol = "amino_acid_change",fillCol="Cohort",txtSize = 4, txtAngle = 90,z=foxo1.dom)


#Supp Fig 1 Plot
foxo1.cleave <- read.table(file="FOXO1_SuppFig1_cleavage_sites.txt",sep='\t',header=T)
foxo1.cleave <- as.data.frame(foxo1.cleave)
foxo1.initiation <- read.table(file="FOXO1_SuppFig1_initiation_sites.txt",sep='\t',header=T)
foxo1.initiation <- as.data.frame(foxo1.initiation)

lolliplot(x=foxo1.cleave,y=foxo1.initiation,txtSize=4,z=foxo1.dom2)


