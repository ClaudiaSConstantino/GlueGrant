

#Joint Modeling Analysis
install.packages("rjags")
install.packages("JM")
install.packages("JMbayes")
library(nlme)
library(survival)
library(JM)
library(rjags)
library(JMbayes)
library(splines)

#Survival sub-model
#I'm not using time independent variables, but i will include age and sex just to obtain a cox model.

coxdata <- subset(clinical.data[,1:3])
coxdata <- cbind(coxdata,ydata_t0)
#change ID names to ID increment values
coxdata <- coxdata %>% mutate(id = cumsum(PATIENT_ID != lag(PATIENT_ID, default="")))
#coxdata$PATIENT_ID <- NULL
#due to an error in data, i will change the time of discharge to the last day in which gene expression were measured
coxdata[44, 4] = 8

#save(coxdata, file = "coxdata.RData")

# Cox model

coxtest <- coxph(Surv(coxdata$time, coxdata$status) ~ SEX + AGE,  data=coxdata, x=TRUE)
summary(coxtest)

coxFinal <- coxph(Surv(coxdata$time, coxdata$status) ~ SEX,  data=coxdata, x=TRUE)
summary(coxFinal)


# JM package

#Rational function fitted
summary(rational.random.no.missing)

jointFitJM_rational <- jointModel(rational.random.no.missing, coxFinal, timeVar = "time", method = "piecewise-PH-GH", lng.in.kn = 5, iter.EM = 200)
summary(jointFitJM_rational)


#Exponential function fitted
summary(exponential.random.no.missing)

jointFitJM_exponential <- jointModel(exponential.random.no.missing, coxFinal, timeVar = "time", method = "piecewise-PH-GH", lng.in.kn = 5, iter.EM = 200)
summary(jointFitJM_exponential)


#Spline function fitted
summary(spline.random.no.missing)

jointFitJM_spline <- jointModel(spline.random.no.missing, coxFinal, timeVar = "time", method = "spline-PH-aGH", lng.in.kn = 5, iter.EM = 200)
summary(jointFitJM_spline)

#Cubic function fitted
summary(cubic.random.no.missing)

jointFitJM_cubic <- jointModel(cubic.random.no.missing, coxFinal, timeVar = "time", method = "piecewise-PH-GH", lng.in.kn = 5, iter.EM = 200)
summary(jointFitJM_cubic)





#JMbayes package


jointFitJMbayes <- jointModelBayes(splineLME, coxFinal, timeVar = "time")

#Rational function fitted
summary(rational.random.no.missing)

jointFitJMbayes_rational <- jointModelBayes(rational.random.no.missing, coxFinal, timeVar = "time")
summary(jointFitJMbayes_rational)


#Exponential function fitted
summary(exponential.random.no.missing)

jointFitJMbayes_exponential <- jointModelBayes(exponential.random.no.missing, coxFinal, timeVar = "time")
summary(jointFitJMbayes_exponential)


#Spline function fitted
summary(spline.random.no.missing)

jointFitJMbayes_spline <- jointModelBayes(spline.random.no.missing, coxFinal, timeVar = "time")
summary(jointFitJMbayes_spline)

#Cubic function fitted
summary(cubic.random.no.missing)

jointFitJMbayes_cubic <- jointModelBayes(cubic.random.no.missing, coxFinal, timeVar = "time")
summary(jointFitJMbayes_cubic)


#Mixed model

#simple function w/ all 12 genes
MixedModelFit <- mvglmer(list(log(NM_013450) ~ time + (time | id),
                              log(AW474434) ~ time + (time | id),
                              log(NM_021730) ~ time + (time | id),
                              log(BG120535) ~ time + (time | id),
                              log(NM_005354) ~ time + (time | id),
                              log(AF279899) ~ time + (time | id),
                              log(BF940270) ~ time + (time | id),
                              log(NM_002600) ~ time + (time | id),
                              log(AW574504) ~ time + (time | id),
                              log(NM_018368) ~ time + (time | id),
                              log(NM_025151) ~ time + (time | id),
                              log(BC000896) ~ time + (time | id)),
                         data = genes.subset2,
                         families = list(gaussian, gaussian,gaussian, gaussian,
                                         gaussian, gaussian, gaussian, gaussian,
                                         gaussian, gaussian, gaussian, gaussian))
summary(MixedModelFit)


CoxFit <- coxph(Surv(coxdata$time, coxdata$status) ~ SEX, data = coxdata, model = TRUE)
summary(CoxFit)

JMFit2 <- mvJointModelBayes(MixedModelFit, CoxFit, timeVar = "time")
summary(JMFit2)




