#Install any of the following packages you don't have (not all are now necessary,
# but at some point they have been, so remove at your own risk)
BiocManager::install(c('epiDisplay',
                       'PoweR',
                       'devtools',
                       'stringi',
                       'survminer',
                       'dplyr',
                       'ggplot2','Rcpp',
                       'purrr',
                       'tidyr',
                       'lmtest',
                       'boot',
                       'colorspace',
                       'coxrt',
                       'pillar',
                       'devEMF',
                       'reshape'))

library('epiDisplay')
library('PoweR')
library('devtools')
library('stringi')
library('survminer')
library('dplyr')
library('ggplot2')
library('Rcpp')
library('purrr')
library('tidyr')
library('lmtest')
library('boot')
library('colorspace')
library('coxrt')
library('pillar')
library('devEMF')
library('reshape')
library('survival')

#Remember if you paste the path in from Windows Explorer to switch the slashes from \ to /
setwd("/path/to/folder")


survdata <- read.csv("*.csv", header = TRUE, row.names=1)

#sBL times (swap in OS.*.3yr or TTP.*.3yr depending on the analysis)
timeOS <- survdata$OS_3yr
eventOS <- survdata$OS_3yr_event

timeTTP <- survdata$TTP_3yr  
eventTTP <- survdata$TTP_3yr_event

#Subset for only cases with data available during the pilot
#survdata <- survdata[survdata$Pilot_37_pBL == 'Yes',]

# survdata$Age.thresh18 <- survdata$Age..if.known.
# survdata$Age.thresh18[survdata$Age..if.known. <= 18] <- c("Young")
# survdata$Age.thresh18[survdata$Age..if.known. > 18] <- c("Old")
# survdata$Age.thresh18 <- as.factor(survdata$Age.thresh18)

#Check the values are now in the survdata object
dim(survdata)
str(survdata) # This tells you what dimensions to put in the next command

#Make sure the subset in the next line is correct for your file. Anything you want to include as a variate
#has to be in the rows included. Cut out the columns you don't want and the survival columns.

#Subset the data to remove descriptor columns - removes an error when printing the univariate results.

ncol(survdata)
covariates <- survdata[,19:ncol(survdata)]
covariates <- colnames(covariates)
covariates


univ_formulasOS <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timeOS, eventOS)~', x)))

univ_formulasTTP <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timeTTP, eventTTP)~', x)))


#Running the Cox model on both OS and TTP datasets - NB: Ensure data isn't mright format
univ_modelsOS <- lapply(univ_formulasOS, function(x){coxph(x, data = survdata)})
univ_modelsTTP <- lapply(univ_formulasTTP, function(x){coxph(x, data = survdata)})


# Extract data
univ_resultsOS <- lapply(univ_modelsOS,
                       function(x){
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (",
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
univ_resultsTTP <- lapply(univ_modelsTTP,
                         function(x){
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR <- paste0(HR, " (",
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, wald.test, p.value)
                           names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                         "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })

# Remember to rename output file here first!

resOS <- t(as.data.frame(univ_resultsOS, check.names = FALSE))
resOS <- as.data.frame(resOS)
resOS

write.table(resOS, "OS_rebute2.txt",sep='\t',quote = F)

resTTP <- t(as.data.frame(univ_resultsTTP, check.names = FALSE))
resTTP <- as.data.frame(resTTP)
resTTP

write.table(resTTP, "TTP_rebute2.txt",sep='\t',quote = F)



#Cohort survival time - OS
surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ 1, data = survdata)
resCOS <- summary(fit, times = c(12))
save.df <- as.data.frame(resCOS[c( "time", "n.risk", "n.event", "surv", "std.err", "lower", "upper")])
save.df.OS

#Cohort survival time - TTP
surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ 1, data = survdata)
resCTTP <- summary(fit, times = c(12))
save.df <- as.data.frame(resCTTP[c( "time", "n.risk", "n.event", "surv", "std.err", "lower", "upper")])
save.df.TTP



#Median OS calulation when survival doesn't go below 50%.
surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$OS.Event.3yr ==1 , data = survdata)
fitOS
plot(fitOS)


#Median TTP calulation when survival doesn't go below 50%.
surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$TTP.Event.3yr ==1 , data = survdata)
fitTTP
plot(fitTTP)


# Multivariate for the RR - forward selection.
TP53_bi_fit <- coxph(Surv(timeTTP, eventTTP) ~
                 survdata$Biallelic_vs_normal +
                 survdata$CNS +
                 survdata$BM + 
                 survdata$LDH_high
               ,data = survdata)
summary(TP53_bi_fit)

TP53_mono_fit <- coxph(Surv(timeTTP, eventTTP) ~
                       survdata$Mono_vs_normal +
                       survdata$CNS +
                       survdata$BM + 
                       survdata$LDH_high
                     ,data = survdata)
summary(TP53_mono_fit)

TP53_any_fit <- coxph(Surv(timeTTP, eventTTP) ~
                         survdata$Any_TP53 +
                         survdata$CNS +
                         survdata$BM + 
                         survdata$LDH_high
                       ,data = survdata)
summary(TP53_any_fit)

TP53_both_fit <- coxph(Surv(timeTTP, eventTTP) ~
                        survdata$Mono_vs_other +
                        survdata$Biallelic_vs_other +
                        survdata$CNS +
                        survdata$BM + 
                        survdata$LDH_high
                      ,data = survdata)
summary(TP53_both_fit)


# eFit2 <- coxph(Surv(time, event) ~
#                  survdata$G_all_chr17.26880986.28944928..Gain +
#                  survdata$GenderF
#                ,data = survdata)
# summary(eFit2)
# lrtest(eFit2,eFit1)
# eFit3 <- coxph(Surv(time, event) ~
#                  survdata$G_all_chr17.26880986.28944928..Gain +
#                  survdata$MCL1.Gain_J19.
#                ,data = survdata)
# summary(eFit3)
# lrtest(eFit3,eFit1)
# eFit4 <- coxph(Surv(time, event) ~
#                  survdata$G_all_chr17.26880986.28944928..Gain +
#                  survdata$MCL1.Gain_J19. +
#                  survdata$GenderF
#                ,data = survdata)
# summary(eFit4)
# lrtest(eFit4,eFit3)
# 
# 
# sFit5 <- coxph(Surv(time, event) ~
#                  survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                  survdata$Age..if.known.
#                ,data = survdata)
# summary(sFit5)
# lrtest(sFit5,sFit1)
# sFit6 <- coxph(Surv(time, event) ~
#                  survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                  survdata$X17p.CNN.LOH +
#                  survdata$X17q.CNN.LOH..manual...10mb.
#                ,data = survdata)
# summary(sFit6)
# lrtest(sFit6,sFit3)
# sFit7 <- coxph(Surv(time, event) ~
#                  survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                  survdata$X17q.CNN.LOH..manual...10mb.+
#                  survdata$No_abn_50kb
#                ,data = survdata)
# summary(sFit7)
# lrtest(sFit7,sFit3)
# sFit8 <- coxph(Surv(time, event) ~
#                  survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                  survdata$X17q.CNN.LOH..manual...10mb. +
#                  survdata$Age..if.known.
#                ,data = survdata)
# summary(sFit8)
# lrtest(sFit8,sFit3):
#   
#   sFit9 <- coxph(Surv(time, event) ~
#                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                    survdata$X17p.CNN.LOH +
#                    survdata$X17q.CNN.LOH..manual...10mb. +
#                    survdata$No_abn_50kb
#                  ,data = survdata)
# summary(sFit9)
# lrtest(sFit9,sFit6)
# 
# sFit10 <- coxph(Surv(time, event) ~
#                   survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                   survdata$X17p.CNN.LOH +
#                   survdata$X17q.CNN.LOH..manual...10mb. +
#                   survdata$Age..if.known.
#                 ,data = survdata)
# summary(sFit10)
# lrtest(sFit10,sFit6)
# 
# sFit11 <- coxph(Surv(time, event) ~
#                   survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
#                   survdata$X17p.CNN.LOH +
#                   survdata$X17q.CNN.LOH..manual...10mb. +
#                   survdata$Age..if.known. +
#                   survdata$No_abn_50kb
#                 ,data = survdata)
# summary(sFit11)
# lrtest(sFit11,sFit10)

#See if Bone Marrow, CSF, CNS, Stage etc add to model
sFitClin1 <- coxph(Surv(time, event) ~
                     survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                     survdata$X17p.CNN.LOH +
                     survdata$X17q.CNN.LOH..manual...10mb. +
                     survdata$Bone.marrow.involvement +
                     survdata$CNS.involvement +
                     survdata$Disease.Stage..St.Jude.
                   ,data = survdata)
summary(sFitClin1)
lrtest(sFitClin1,sFit10)

lrtest(sFit2,sFit1)
lrtest(sFit3,sFit1)
lrtest(sFit4,sFit1)
lrtest(sFit5,sFit1)
lrtest(sFit6,sFit2)

#Multivariate for the RR - backward selection.

bFit1 <- coxph(Surv(time, event) ~
                 survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                 survdata$X17p.CNN.LOH +
                 survdata$X17q.CNN.LOH..manual...10mb. +
                 survdata$Age..if.known. +
                 survdata$No_abn_50kb
               ,data = survdata)
summary(bFit1)
#Removing one factor each time
bFit2 <- coxph(Surv(time, event) ~
                 survdata$X17p.CNN.LOH +
                 survdata$X17q.CNN.LOH..manual...10mb. +
                 survdata$Age..if.known. +
                 survdata$No_abn_50kb
               ,data = survdata)
summary(bFit2)
lrtest(bFit2,bFit1)

bFit3 <- coxph(Surv(time, event) ~
                 survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                 survdata$X17q.CNN.LOH..manual...10mb. +
                 survdata$Age..if.known. +
                 survdata$No_abn_50kb
               ,data = survdata)
summary(bFit3)
lrtest(bFit3,bFit1)

bFit4 <- coxph(Surv(time, event) ~
                 survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                 survdata$X17p.CNN.LOH +
                 survdata$Age..if.known. +
                 survdata$No_abn_50kb
               ,data = survdata)
summary(bFit4)
lrtest(bFit4,bFit1)

bFit5 <- coxph(Surv(time, event) ~
                 survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                 survdata$X17p.CNN.LOH +
                 survdata$X17q.CNN.LOH..manual...10mb. +
                 survdata$No_abn_50kb
               ,data = survdata)
summary(bFit5)
lrtest(bFit5,bFit1)

bFit6 <- coxph(Surv(time, event) ~
                 survdata$v9_GISTIC2_chr3.195718766.198022430..Gain + 
                 survdata$X17p.CNN.LOH + 
                 survdata$X17q.CNN.LOH..manual...10mb. + 
                 survdata$Age..if.known. 
               ,data = survdata)
summary(bFit6)
lrtest(bFit6,bFit1)


summary(sFit6)
HR <- round(exp(coef(sFit)), 2)
CI <- round(exp(confint(sFit)), 2)
P <- round(coef(summary(sFit))[,5], 3)
# Names the columns of CI
colnames(CI) <- c("Lower", "Higher")
# Bind columns together as dataset
table2 <- as.data.frame(cbind(HR, CI, P))
table2
#  #   #write.csv(table2, "RR_MultivariateForTier3_analysis2.csv")

######## BIALLELIC STUFF ########
TP53Fit0 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain
                  ,data = survdata)
summary(TP53Fit0)

TP53Fit1 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$X17p.CNN.LOH
                  ,data = survdata)
summary(TP53Fit1)
lrtest(TP53Fit1,TP53Fit0)

TP53Fit2 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$X17q.CNN.LOH..manual...10mb.
                  ,data = survdata)
summary(TP53Fit2)
lrtest(TP53Fit2,TP53Fit0)

TP53Fit3 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$TP53.Biallelic
                  ,data = survdata)
summary(TP53Fit3)
lrtest(TP53Fit3,TP53Fit0)

TP53Fit4 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$TP53.Mutation..Y.N.
                  ,data = survdata)
summary(TP53Fit4)
lrtest(TP53Fit4,TP53Fit0)

TP53Fit5 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$X17p.CNN.LOH +
                    survdata$X17q.CNN.LOH..manual...10mb.
                  ,data = survdata)
summary(TP53Fit5)
lrtest(TP53Fit5,TP53Fit1)

TP53Fit6 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$X17p.CNN.LOH +
                    survdata$TP53.Biallelic
                  ,data = survdata)
summary(TP53Fit6)
lrtest(TP53Fit6,TP53Fit1)

TP53Fit7 <- coxph(Surv(time, event) ~
                    survdata$v9_GISTIC2_chr3.195718766.198022430..Gain +
                    survdata$X17p.CNN.LOH +
                    survdata$TP53.Mutation..Y.N.
                  ,data = survdata)
summary(TP53Fit7)
lrtest(TP53Fit7,TP53Fit1)

#Survival curves

#Cohort Plot
surv <- Surv(time,event)
fit <- survfit(surv ~ 1, data = survdata)
names(fit$strata)
summary(fit)
plot(fit)
pdf ("OS_Cohort_Malawi_all.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "Cohort.",
           legend.labs = c("Cohort"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (years)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           #pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = T,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

##########################################################################################

#Any TP53 TTP

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$Any_TP53, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_Any_TP53.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Abnormality",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (years)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

#17p CNN-LOH TTP

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$X17p_CNN.LOH, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_17p_CNN-LOH.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "17p CNN-LOH",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

# 17p Loss

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$X17p_Loss, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_17p_Loss.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "17p Loss",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

#TP53 Mutation TTP

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$Mutation, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_TP53_mutation.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Mutation",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()


#Biallelic vs normal TTP

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$Biallelic_vs_normal, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_Biallelic_norm.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Biallelic",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

#Monoallelic vs others TTP

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$Monoallelic_vs_normal, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_Monoallelic_norm.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Monoallelic",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

# TP53 Split Plot TTP

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$Biallelic.Monoallelic.Normal_Exc_Silent, data = survdata)
names(fit$strata)

plot(fit)
pdf ("TTP_TP53_Split.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Status",
           legend.labs = c("Biallelic","Monoallelic","Normal"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Time to Progression Survival %",
           # Add p-value and tervals
           pval = F,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#CF3A3A", "#ffa501ff","black"),
           linetype = c("solid","solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

############################################ OS ###

#Any TP53 OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$Any_TP53, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_Any_TP53.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Abnormality",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (years)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

#17p CNN-LOH OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$X17p_CNN.LOH, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_17p_CNN-LOH.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "17p CNN-LOH",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

# 17p Loss OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$X17p_Loss, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_17p_Loss.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "17p Loss",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

#TP53 Mutation OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$Mutation, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_TP53_mutation.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Mutation",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()


#Biallelic vs normal OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$Biallelic_vs_normal, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_Biallelic_norm.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Biallelic",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

#Monoallelic vs others OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$Monoallelic_vs_normal, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_Monoallelic_norm.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Monoallelic",
           legend.labs = c("No","Yes"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

# TP53 Split Plot OS

surv <- Surv(timeOS,eventOS)
fit <- survfit(surv ~ survdata$Biallelic.Monoallelic.Normal_Exc_Silent, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_TP53_Split.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "TP53 Status",
           legend.labs = c("Biallelic","Monoallelic","Normal"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = F,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#CF3A3A", "#ffa501ff","black"),
           linetype = c("solid","solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()





#Cohort Plot
surv <- Surv(time,event)
fit <- survfit(surv ~ 1, data = survdata)
names(fit$strata)
summary(fit)
plot(fit)
pdf ("OS_Cohort_Malawi_all.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "Cohort.",
           legend.labs = c("Cohort"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (years)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           #pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = T,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()


#Plot Gender

surv <- Surv(time,event)
fit <- survfit(surv ~ survdata$Gender, data = survdata)
names(fit$strata)

plot(fit)
pdf ("OS_Gender.pdf")
#ggsurvplot(fit,font.legend = c(, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
ggsurvplot(fit,
           # Change legends: title & labels
           font.x = c(18,"bold","black"),
           font.y = c(18,"bold","black"),
           font.legend = c(14,"bold","black"),
           legend.title = "Gender",
           legend.labs = c("Female","Male"),
           font.tickslab = c("14","plain","black"),
           xlab="Follow up (months)",
           ylab="Overall Survival %",
           # Add p-value and tervals
           pval = TRUE,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("black", "#CF3A3A"),
           linetype = c("solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()


#For 3 options -RR





#Log Rank Test
survdiff(Surv(time,event) ~ X17p.CNN.LOH,data=survdata)

#Can't make the above command turn the risk table labels black...

# Section below is for the old favourite of dashed and solid black and white lines in KMs.
# surv <- Surv(time,event)
# fit <- survfit(surv ~ survdata$v9_GISTIC2_chr3.195718766.198022430..Gain, data = survdata)
# names(fit$strata)
# #names(fit$strata) <- c("No 17p","Yes 17p")
# plot(fit)
# pdf ("RR_3q29_gain_all.file")
# #ggsurvplot(fit,font.legend = c(6, "plain", "black"),linetype = c("solid","dashed"),palette = "grey",risk.table = "absolute")
# ggsurvplot(fit,
#            # Change legends: title & labels
#            legend.title = "TP53 Biallelic",
#            legend.labs = c("No", "Yes"),
#            xlab="Overall Survival in years",
#            # Add p-value and tervals
#            pval = FALSE,
#            
#            conf.int = FALSE,
#            # Add risk table
#            risk.table = TRUE,
#            tables.height = 0.2,
#            tables.theme = clean_theme(),
#            
#            # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
#            # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
#            palette = c("#0f0e0c", "#0f0e0c"),
#            linetype = c("solid","dashed"),
#            ggtheme = theme_survminer() # Change ggplot2 theme
# )
# dev.off()



surv <- Surv(time,event)
# fit <- survfit(surv ~ 1 , data = survdata)
# res <- summary(fit, times = c(1,2,3))
# save.df <- as.data.frame(res[c( "time", "n.risk", "n.event", "surv", "std.err", "lower", "upper")])
# save.df
# write.csv(save.df, "RR_tier3_KmEstimates.csv")
