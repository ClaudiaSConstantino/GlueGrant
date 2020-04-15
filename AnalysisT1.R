# BY CLÁUDIA CONSTANTINO

#Analysis of gene data at the first day post injury (162 patients included)


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
data_t1 <- patient.data[(patient.data$SAMPLE_STUDYSTART_DAYS == "1"),]
data_t1$ICU_VENTDAYS <- NULL
data_t1$SAMPLE_STUDYSTART_DAYS <- NULL

#select the common patients, the ones that have microarray data at time = 1
common <- intersect(data_t1$PATIENT_ID, clinical.data$PATIENT_ID)  
clinical.data.t1 <- clinical.data[clinical.data$PATIENT_ID %in% common,]


data_t1$DEATH_OR_DISCHARGE_DAY_SNC_INJ <- clinical.data.t1$DEATH_OR_DISCHARGE_DAY_SNC_INJ
data_t1$EVENT <- clinical.data.t1$DISCHARGE_DAY_SNC_INJ

#create vector with 1 and 0 corresponding to event
data_t1$EVENT <- ifelse(data_t1$EVENT > 1 , 1, data_t1$EVENT)
data_t1$EVENT[is.na(data_t1$EVENT)] <- 0 #substitue NA's for 0

ydata_t1 <- cbind(data_t1$DEATH_OR_DISCHARGE_DAY_SNC_INJ, data_t1$EVENT)
colnames(ydata_t1) <- c("time", "status")

#to see in which index is the follwoing varibles
grep("DEATH_OR_DISCHARGE_DAY_SNC_INJ", colnames(data_t1))
grep("EVENT", colnames(data_t1))
data_t1[1:10, 1:10]

genes_t1 <- as.data.frame(data_t1[c(9:54683)])
genes_t1 <- as.matrix(genes_t1)


#PACKAGE: GLMNET

ydata <- as.matrix(ydata_t1)

cvfit = cv.glmnet(genes_t1, ydata_t1, family = "cox")

fit <- glmnet(genes_t1, ydata_t1, 
              alpha = 0.8,
              lambda = cvfit$lambda.min,
              family = "cox")

elastic.cox.t1 <- coef(fit, s = 'lambda.1se')[,1] %>% {.[.!=0]} #to select just variables with coef different from zero
length(elastic.cox.t1)
names(elastic.cox.t1)
genes.selected.t1 <- names(elastic.cox.t1)







