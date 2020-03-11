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

setwd("C:/Users/U0034370/OneDrive - Teesside University/Paper_Rebuttal/Survival")
full_input <- read.csv("Updated_Survival_Master_19022020_v1.csv", header = TRUE, row.names=1)

#Create a directory for your files

dir.create('Include_all_pvalue_without_conf_int')

setwd("C:/Users/U0034370/OneDrive - Teesside University/Paper_Rebuttal/Survival/Include_all_pvalue_without_conf_int")

### Set your cohort

#All
survdata <- full_input
dim(survdata)

#Include_90 all with clinical data
survdata <- full_input[full_input$Include_90 == 'Yes',]
dim(survdata)

#Include_81 all with clinical data and FAB/LMB96 treated
survdata <- full_input[full_input$Include_81 == 'Yes',]
dim(survdata)

#Include_BL with clinical data n=
survdata <- full_input[full_input$Diagnosis == 'BL',]
dim(survdata)

#Include_BL with clinical data and FAB/LMB96 treated n=
survdata <- full_input[full_input$Diagnosis == 'BL' & full_input$Include_81 == 'Yes',]
dim(survdata)

#sBL times (swap in OS.*.3yr or TTP.*.3yr depending on the analysis)
timeOS <- survdata$OS_3yr
eventOS <- survdata$OS_3yr_event

timeTTP <- survdata$TTP_3yr  
eventTTP <- survdata$TTP_3yr_event

### Check the data is showing what you want ###

dim(survdata)
str(survdata)


### Set up what format you want for the KM plots

pval <- TRUE
conf.int <- F


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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#CF3A3A", "#ffa501ff","#548235"),
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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
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
           pval = pval,
           risk.table.col = "black",
           fontsize = 6,
           
           conf.int = conf.int,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = clean_theme(),
           
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#CF3A3A", "#ffa501ff","#548235"),
           linetype = c("solid","solid","solid"),
           ggtheme = theme_survminer() # Change ggplot2 theme
)
dev.off()

# survdiff(Surv(timeTTP,eventTTP) ~ Any_TP53,data=survdata)
# 
# data.survdiff <- survdiff(Surv(timeTTP,eventTTP) ~ Any_TP53,data=survdata)
# p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
# HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
# up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
# low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
# HR
# low95
# up95

surv <- Surv(timeTTP,eventTTP)
fit <- survfit(surv ~ survdata$Biallelic_vs_normal, data = survdata)
fit
