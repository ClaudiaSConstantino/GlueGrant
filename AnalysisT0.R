if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("glmSparseNet")

library(dplyr)
library(ggplot2)
library(survival)
library(loose.rock)
library(futile.logger)
#
library(glmSparseNet)
library(glmnet)
#
# Some general options for futile.logger the debugging package
.Last.value <- flog.layout(layout.format('[~l] ~m'))
.Last.value <- loose.rock::show.message(FALSE)
# Setting ggplot2 default theme as minimal
theme_set(ggplot2::theme_minimal())


library(tidyverse)
library(readxl)

#Load file "GluegrantData.RData"

#Load clinical data
clinical.data <- read_excel(
  "TRDB_TRAUMA-PATIENT_CLIN_DATA_RPT_20090126_195655.xls")

clinical.data <- subset(clinical.data, select = c(PATIENT_ID, AGE, SEX, RACE, HEIGHT_CM, WEIGHT, 
                                                  ICU_DAYS, DISCHARGE_DAY_SNC_INJ, DEATH_OR_DISCHARGE_DAY_SNC_INJ))

head(clinical.data)
head(patient.data)

#write.table(patient.data, "patient.data.txt")


#select in data t0 for a simple cox analysis
data_t0 <- patient.data[(patient.data$SAMPLE_STUDYSTART_DAYS == "0"),]
data_t0$ICU_VENTDAYS <- NULL
data_t0$SAMPLE_STUDYSTART_DAYS <- NULL
data_t0$DEATH_OR_DISCHARGE_DAY_SNC_INJ <- clinical.data$DEATH_OR_DISCHARGE_DAY_SNC_INJ
data_t0$EVENT <- clinical.data$DISCHARGE_DAY_SNC_INJ

data_t0$EVENT <- ifelse(data_t0$EVENT > 1 , 1, data_t0$EVENT)
data_t0$EVENT[is.na(data_t0$EVENT)] <- 0 #substitue NA's for 0

ydata_t0 <- cbind(data_t0$DEATH_OR_DISCHARGE_DAY_SNC_INJ, data_t0$EVENT)
colnames(ydata_t0) <- c("time", "status")


genes_t0 <- as.data.frame(data_t0[c(9:54683)])
genes_t0 <- as.matrix(genes_t0)



#TO COMPARE SELECTED GENES WITH OTHER TIMES (T1, T3...)

ydata_t0 <- as.matrix(ydata_t0)

cvfit = cv.glmnet(genes_t0, ydata_t0, family = "cox")

fit <- glmnet(genes_t0, ydata_t0, 
              alpha = 0.8,
              lambda = cvfit$lambda.min,
              family = "cox")

elastic.cox.t0 <- coef(fit, s = 'lambda.1se')[,1] %>% {.[.!=0]} #to select just variables with coef different from zero
length(elastic.cox.t0)
names(elastic.cox.t0)
genes.selected.t0 <- names(elastic.cox.t0)


#NOT INCLUDED YET
#TESTS WITH DATA SPLIT INTO TEST/TRAIN

#split 70% for train and 30% for test 
index <- sample(nrow(genes_t0), 118)
datatreino <- as.matrix(genes_t0[index, ])
datateste <- as.matrix(genes_t0[-index, ])
vectortreino <- ydata[index, ]
vectorteste <- ydata[-index, ]




#PACKAGE: glmSparseNet
fitted <- cv.glmHub(datatreino, Surv(vectortreino$time, vectortreino$status),
                    family  = 'cox',
                    lambda = buildLambda(1),
                    network = 'correlation', 
                    network.options = networkOptions(cutoff = .6, 
                                      min.degree = .2))

summary(fitted)

coefs.v <- coef(fitted, s = 'lambda.min')[,1] %>% { .[. != 0]}

# to give the genes selected and the respective coeficients
coefs.v %>% { 
  data.frame(gene.name   = names(.),
             coefficient = .,
             stringsAsFactors = FALSE)
} %>%
  arrange(gene.name) %>%
  knitr::kable()

summary(coefs.v)
length(coefs.v)
names(coefs.v)

#to select in the data the genes given by cv.glmHub method 
coef.names <- as.vector(names(coefs.v))
#genes_t0 <- as.data.frame(genes_t0)
names.use <- names(genes_t0)[(names(genes_t0) %in% coef.names)]
genes_t0_subset <- genes_t0[,names.use]



#to compute survival curve
#TREINO
graphtrain <- separate2GroupsCox(as.vector(coefs.v), 
                   datatreino[, names(coefs.v)], 
                   vectortreino, 
                   plot.title = 'Full dataset', legend.outside = FALSE)


#Test
separate2GroupsCox(as.vector(coefs.v), 
                   datateste[, names(coefs.v)], 
                   vectorteste, 
                   plot.title = 'Full dataset', legend.outside = FALSE)




#univariate cox analysis  
cox.NM <- coxph(Surv(ydata$time, ydata$status) ~ NM_006276, data = genes_t0_subset)
summary(cox.NM)







#PACKAGE: GLMNET


vectortreino <- as.matrix(vectortreino)
cvfit = cv.glmnet(datatreino, vectortreino, family = "cox")

fit <- glmnet(datatreino, vectortreino, 
              alpha = 1,
              lambda = cvfit$lambda.min,
              family = "cox")

elastic.cox <- coef(fit, s = 'lambda.1se')[,1] %>% {.[.!=0]} #to select just variables with coef different from zero
length(elastic.cox)
names(elastic.cox)


#to compute survival curve
#TREINO
vectortreino <- as.data.frame(vectortreino)
separate2GroupsCox(as.vector(elastic.cox), 
                                 datatreino[, names(elastic.cox)], 
                                 vectortreino, 
                                 plot.title = 'Full dataset', legend.outside = FALSE)


#TEST
separate2GroupsCox(as.vector(elastic.cox), 
                   datateste[, names(elastic.cox)], 
                   vectorteste, 
                   plot.title = 'Full dataset', legend.outside = FALSE)

#multivariate cox analysis  

#to select in the data the genes given by cv.glmHub method 
coef.names <- as.vector(names(elastic.cox))
#genes_t0 <- as.data.frame(genes_t0)
names.use <- names(genes_t0)[(names(genes_t0) %in% coef.names)]
genes_t0_subset <- genes_t0[,names.use]


multicox <- coxph(Surv(ydata$time, ydata$status) ~  NM_002379 + AF151074 + N37081 +
                  AI022882 +  NM_024873 + BC005131  + N63748 + BE675685 + AA808444 + AW014730 + 
                  AI076012 +  AW086493 + AK095151 + BG120535 + AK098018, data = genes_t0)

summary(multicox)




