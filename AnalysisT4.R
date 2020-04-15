# BY CLÁUDIA CONSTANTINO

#Analysis of gene data at the fourth day post injury (111 patients included)

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


#select in data t0 for a simple cox analysis
data_t4 <- patient.data[(patient.data$SAMPLE_STUDYSTART_DAYS == "4"),]
data_t4$ICU_VENTDAYS <- NULL
data_t4$SAMPLE_STUDYSTART_DAYS <- NULL

#select the common patients, the ones that have microarray data at time = 4
common <- intersect(data_t4$PATIENT_ID, clinical.data$PATIENT_ID)  
clinical.data.t4 <- clinical.data[clinical.data$PATIENT_ID %in% common,]


data_t4$DEATH_OR_DISCHARGE_DAY_SNC_INJ <- clinical.data.t4$DEATH_OR_DISCHARGE_DAY_SNC_INJ
data_t4$EVENT <- clinical.data.t4$DISCHARGE_DAY_SNC_INJ




#create vector with 1 and 0 corresponding to event
data_t4$EVENT <- ifelse(data_t4$EVENT > 1 , 1, data_t4$EVENT)
data_t4$EVENT[is.na(data_t4$EVENT)] <- 0 #substitue NA's for 0
data_t4$EVENT

ydata_t4 <- cbind(data_t4$DEATH_OR_DISCHARGE_DAY_SNC_INJ, data_t4$EVENT)
colnames(ydata_t4) <- c("time", "status")

#to see in which index is the follwoing varibles
grep("DEATH_OR_DISCHARGE_DAY_SNC_INJ", colnames(data_t4))
grep("EVENT", colnames(data_t4))
data_t4[1:10, 1:10]

genes_t4 <- as.data.frame(data_t4[c(9:54683)])
genes_t4 <- as.matrix(genes_t4)


#PACKAGE: GLMNET

ydata_t4 <- as.matrix(ydata_t4)

cvfit = cv.glmnet(genes_t4, ydata_t4, family = "cox")

fit <- glmnet(genes_t4, ydata_t4, 
              alpha = 0.8,
              lambda = cvfit$lambda.min,
              family = "cox")

elastic.cox.t4 <- coef(fit, s = 'lambda.1se')[,1] %>% {.[.!=0]} #to select just variables with coef different from zero
length(elastic.cox.t4)
names(elastic.cox.t4)
genes.selected.t4 <- names(elastic.cox.t4)









