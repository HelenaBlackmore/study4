#####################################################################
### 3. Model-based causal mediation analysis
#####################################################################
### 3.1. Estimation of the average causal mediation effects
library("mediation")
set.seed(2014)
data("framing", package = "mediation")

med.fit <- lm(emo ~ treat + age + educ + gender + income, data = framing)
out.fit <- glm(cong_mesg ~ emo + treat + age + educ + gender + income,
  data = framing, family = binomial("probit"))

med.out <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo",
  robustSE = TRUE, sims = 100)
summary(med.out)

med.out <- mediate(med.fit, out.fit, boot = TRUE, treat = "treat",
  mediator = "emo", sims = 100)
summary(med.out)

med.fit <- lm(emo ~ treat + age + educ + gender + income, data = framing)
out.fit <- glm(cong_mesg ~ emo * treat + age + educ + gender + income,
  data = framing, family = binomial("probit"))
med.out <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo",
  robustSE = TRUE, sims = 100)
summary(med.out)

test.TMint(med.out, conf.level = 0.95)


### 3.2. Moderated mediation
med.fit <- lm(emo ~ treat * age + educ + gender + income, data = framing)
out.fit <- glm(cong_mesg ~ emo + treat * age + emo * age + educ + gender +
  income, data = framing, family = binomial("probit"))

med.age20 <- mediate(med.fit, out.fit, treat = "treat",
  mediator = "emo", covariates = list(age = 20), sims = 100)
med.age60 <- mediate(med.fit, out.fit, treat = "treat",
  mediator = "emo", covariates = list(age = 60), sims = 100)

summary(med.age20)
summary(med.age60)

med.init <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo", 
  sims = 2)
test.modmed(med.init, covariates.1 = list(age = 20),
  covariates.2 = list(age = 60), sims = 100)

### 3.3. Non-binary treatment variables
med.fit <- lm(emo ~ cond + age + educ + gender + income, data = framing)
out.fit <- glm(cong_mesg ~ emo + cond + age + educ + gender + income,
  data = framing, family = binomial("probit"))
med23.out <- mediate(med.fit, out.fit, treat = "cond", mediator = "emo",
  control.value = 2, treat.value = 3, sims = 100)
summary(med23.out)

med14.out <- mediate(med.fit, out.fit, treat = "cond", mediator = "emo",
  control.value = 1, treat.value = 4, sims = 100)
summary(med14.out)

### 3.4. Sensitivity analysis for sequential ignorability
med.fit <- lm(emo ~ treat + age + educ + gender + income, data = framing)
out.fit <- glm(cong_mesg ~ emo + treat + age + educ + gender + income,
  data = framing, family = binomial("probit"))
med.out <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo",
  robustSE = TRUE, sims = 100)
sens.out <- medsens(med.out, rho.by = 0.1, effect.type = "indirect", sims = 100)
summary(sens.out)

plot(sens.out, sens.par = "rho", main = "Anxiety", ylim = c(-0.2, 0.2))
plot(sens.out, sens.par = "R2", r.type = "total", sign.prod = "positive")

#####################################################################
### 4. Causal mediation analysis of multilevel data
#####################################################################
### 4.1. Group-level treatment and individual-level mediator
data("student", package = "mediation")
library("lme4")
med.fit <- glmer(attachment ~ catholic + gender + income + pared + (1 | SCH_ID),
  family = binomial(link = "logit"), data = student)
out.fit <- glmer(fight ~ catholic*attachment +
  gender + income + pared + (1 + attachment | SCH_ID),
  family = binomial(link = "logit"), data = student)

med.out <- mediate(med.fit, out.fit, treat = "catholic", mediator = "attachment",
  sims = 100)
summary(med.out)

### 4.2. Group-level treatment and mediator
data("school", package = "mediation")
med.fit <- lm(smorale ~  free + catholic + coed, data = school)
out.fit <- lmer(late ~ free + smorale + catholic + coed +
  gender + income + pared + (1 | SCH_ID), data = student)

med.out <- mediate(med.fit, out.fit, treat = "free", mediator = "smorale",
  control.value = 3, treat.value = 4, sims = 100)
summary(med.out)

#####################################################################
### 5. Design-based causal mediation analysis
#####################################################################
### 5.1. Single experiment design
framing$english <- as.numeric(framing$english)
framing$anx <- as.numeric(framing$anx)
sed.est <- mediate.sed("english", "anx", "treat", data = framing, SI = TRUE,
  boot = TRUE, sims = 100)
summary(sed.est)

### 5.2. Parallel design
data("boundsdata", package = "mediation")
pd <- mediate.pd("out", "med", "ttt", "manip", boundsdata,
  NINT = TRUE, sims = 100, conf.level = 0.95)
summary(pd)
pd1 <- mediate.pd("out", "med", "ttt", "manip", boundsdata,
                 NINT = FALSE)
summary(pd1)

### 5.3. Parallel encouragement design
data("boundsdata", package = "mediation")
ped <- mediate.ped("out.enc", "med.enc", "ttt", "enc", boundsdata)
summary(ped)

### 5.4. Crossover encouragement design
data("CEDdata", package = "mediation")
ced <- mediate.ced("Y2", "M1", "M2", "T1", "Z", CEDdata, sims = 100)
summary(ced)


#####################################################################
### 6. Analysis of causally dependent multiple mechanisms
#####################################################################
### 6.2. Single experiment design
Xnames <- c("age", "educ", "gender", "income")
m.med <- multimed(outcome = "immigr", med.main = "emo", med.alt = "p_harm",
  treat = "treat", covariates = Xnames, data = framing, sims = 100)
summary(m.med)

plot(m.med, type = "point")
plot(m.med, type = c("sigma", "R2-total"), tgroup = c("treated", "control"))

### 6.3. Parallel design
m.med.para <- multimed(outcome = "out", med.main = "med", treat = "ttt",
  experiment = "manip", design = "parallel", data = boundsdata, sims = 100)
summary(m.med.para)

#####################################################################
### 7. Causal mediation analysis with treatment nNoncompliance
#####################################################################
### 7.2. Illustration
data("jobs", package = "mediation")
a <- lm(comply ~ treat + sex + age + marital + nonwhite + educ + income, 
  data = jobs)
b <- glm(job_dich ~ comply + treat + sex + age + marital + 
  nonwhite + educ + income, data = jobs, family = binomial)
c <- lm(depress2 ~ job_dich * (comply + treat) + sex + age + marital + 
  nonwhite + educ + income, data = jobs)

out <- ivmediate(a, b, c, sims = 100, boot = FALSE, 
  enc = "treat", treat = "comply", mediator = "job_dich")
summary(out)
plot(out) 
