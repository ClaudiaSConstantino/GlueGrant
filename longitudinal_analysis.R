
#LME and NLME with the 12genes subset "genes.data.subset2"

library(nlme)
library(splines)
library(dplyr)
library(ggplot2)

#linear mixed-effects models (LME)

genes.subset2 <- genes.data.subset2 %>%
  select(PATIENT_ID, MICROARRAY_TIME, everything())

colnames(genes.subset2)[1] <- "patient"
colnames(genes.subset2)[2] <- "time"

#change ID names to ID increment values
genes.subset2 <- genes.subset2 %>% mutate(id = cumsum(patient != lag(patient, default="")))

#save(genes.subset2, file = "genes.subset2.RData")
#save(genes.data.subset2, file = "genes.data.subset2.RData")


# LME Random intercept model
lmeFit.interc <- lme( NM_013450 ~ time, random = ~ 1 | id , data=genes.subset2) 
summary(lmeFit.interc)

# LME Random intercept and random slope model
lmeFit.slope <- lme( NM_013450 ~ time, random = ~ time | id , data=genes.subset2) 
summary(lmeFit.slope)

#deal with the missing data is to omit the missing values, 
#this can be done since the LME model does not require for 
#the measurements to be equally spaced

# RATIONAL function formula 
control = nlmeControl(pnlsTol = x, msVerbose = TRUE)
# Obtain optimal tuning parameters using the *nlme* routine
rational.nlme <- nlme(model = log(NM_013450) ~ 1/(time + a)^g + i,
                      fixed =  a + g + i ~ 1,
                      random = i ~ 1 | id,
                      start = c(a = 1, g = 1, i = 1),
                      control = control,
                      data = genes.subset2)

# Compose the fixed and random effects formulas
fixed.rational.omit <- as.formula(sprintf("log(NM_013450) ~ I(1/(time + %f)^%f) + 1",
                                          coef(rational.nlme)[1,"a"], 
                                          coef(rational.nlme)[1,"g"]))

random.rational.omit <- as.formula(sprintf("~ I(1/(time + %f)^%f) + 1 | id",
                                           coef(rational.nlme)[1,"a"], 
                                           coef(rational.nlme)[1,"g"]))

rational.random.no.missing <- lme(fixed = fixed.rational.omit, 
                                  random = random.rational.omit,
                                  method = "ML", data = genes.subset2,
                                  control = lmeControl(opt = "optim"))

summary(rational.random.no.missing)


#LOOP
colNames <- names(genes.subset2)[3:14]
control = nlmeControl(pnlsTol = x, msVerbose = TRUE)
teste = list()
for(x in colNames){
  n <- noquote(x)
  # Obtain optimal tuning parameters using the *nlme* routine
  rational.nlme <- nlme(model = log(n) ~ 1/(time + a)^g + i,
                        fixed =  a + g + i ~ 1,
                        random = i ~ 1 | id,
                        start = c(a = 1, g = 1, i = 1),
                        control = control,
                        data = genes.subset2)
  
  # Compose the fixed and random effects formulas
  fixed.rational.omit <- as.formula(sprintf("log(gene) ~ I(1/(time + %f)^%f) + 1",
                                            coef(rational.nlme)[1,"a"], 
                                            coef(rational.nlme)[1,"g"]))
  
  random.rational.omit <- as.formula(sprintf("~ I(1/(time + %f)^%f) + 1 | id",
                                             coef(rational.nlme)[1,"a"], 
                                             coef(rational.nlme)[1,"g"]))
  
  rational.random.no.missing <- lme(fixed = fixed.rational.omit, 
                                    random = random.rational.omit,
                                    method = "ML", data = genes.subset2,
                                    control = lmeControl(opt = "optim"))
  teste[[i]] = rational.random.no.missing

}


# EXPONENCIAL function formula
exponential.nlme <- nlme(NM_013450 ~ exp(-time*a) + int,
                                fixed = a + int ~ 1,
                                random = int ~ 1 | id,
                                data = genes.subset2,
                                control = control,
                                start = c( a = 1, int = 3))

fixed.exponential.omit <- as.formula(sprintf("log(NM_013450) ~ I(exp(-time*%f)) + 1",
                                             coef(exponential.nlme)$a[1]))

random.exponential.omit <- as.formula(sprintf("~ I(exp(-time*%f)) + 1 | id",
                                              coef(exponential.nlme)$a[1]))

exponential.random.no.missing <- lme(fixed = fixed.exponential.omit,
                                     random = random.exponential.omit,
                                     method = "ML", data = genes.subset2,
                                     control = lmeControl(opt = "optim"))


summary(exponential.random.no.missing)


# SPLINE function formula
spline.random.no.missing <- lme(fixed = log(NM_013450) ~ ns(time,3) + 1,
                                random = ~ ns(time, 3) + 1 | id,
                                method = "ML", data = genes.subset2,
                                control = lmeControl(opt = "optim"))


summary(spline.random.no.missing)


# CUBIC function formula

# Note that all the random effects except b0 were removed from the 
#Cubic model due to convergence issues.

cubic.random.no.missing1 <- lme(fixed = log(NM_013450) ~ poly(time,3) + 1,
                               random = ~ 1 | id,
                               method = "ML", data = genes.subset2,
                               control = lmeControl(opt = "optim"))
cubic.random.no.missing2 <- lme(fixed = log(AW474434) ~ poly(time,3) + 1,
                               random = ~ 1 | id,
                               method = "ML", data = genes.subset2,
                               control = lmeControl(opt = "optim"))
cubic.random.no.missing3 <- lme(fixed = log(NM_021730) ~ poly(time,3) + 1,
                                random = ~ 1 | id,
                                method = "ML", data = genes.subset2,
                                control = lmeControl(opt = "optim"))

summary(cubic.random.no.missing1)
summary(cubic.random.no.missing2)
summary(cubic.random.no.missing3)




## Plot the fits


data.plot <- genes.subset2
#data.plot$lmeFit.interc <- predict(lmeFit.interc)
#data.plot$lmeFit.slope <- predict(lmeFit.slope)
data.plot$spline.fit <- predict(spline.random.no.missing)
data.plot$exponential.random.no.missing <- predict(exponential.random.no.missing)
data.plot$rational.random.no.missing <- predict(rational.random.no.missing)
data.plot$cubic.random.no.missing <- predict(cubic.random.no.missing)

#data.plot[data.plot$Death.Status == 0 , "Death.Status"] <- "Alive"
#data.plot[data.plot$Death.Status == 1 , "Death.Status"] <- "Dead"



lmegraphs <- ggplot(data = data.plot, aes(x = time, y = log(NM_013450), group = id)) + 
              geom_point(color="black", size=0.3) +
              geom_line(color="black", size=0.5) + 
              #geom_line(aes(y = lmeFit.interc, color = "a")) + 
              #geom_line(aes(y = lmeFit.slope, color = "b")) + 
              geom_line(aes(y = cubic.random.no.missing, color = "c")) + 
              geom_line(aes(y = rational.random.no.missing, color = "d")) + 
              geom_line(aes(y = exponential.random.no.missing, color = "e")) + 
              geom_line(aes(y = spline.fit, color = "f")) + 
              theme_bw()+
              scale_color_manual(name = "Model Fit", 
                                 labels = c("Cubic Model", "Rational Model", "Exponential Model",
                                            "Spline Model"),
                                 values = c("c" = "magenta3", "d" = "deepskyblue", "e" = "red2",
                                            "f" = "green3")) + 
              facet_wrap(~id) +
              scale_x_continuous(breaks = c(0,10,20,30)) + 
              theme(legend.position = "bottom")

print(lmegraphs)
ggsave("lmegraphs.png", width = 10, height = 15)
########








