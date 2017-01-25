library(lme4)
library(lsmeans)
library(data.table)
library(nlme)

# load data
library(readr)

scan_behavior = read_delim("~/Box Sync/Project Trust Game/data/processed/scan_behavior.txt",
                            "\t", escape_double = FALSE, trim_ws = TRUE)
#beha_behavior = read_delim("~/Box Sync/Project Trust Game/data/processed/beha_behavior.txt",
#                           "\t", escape_double = FALSE, trim_ws = TRUE)
#b = data.table(beha_behavior)
b = data.table(scan_behavior)
#attach(b)
View(b)

#lag variable: previous trustee decision
b[, pt_decision:=c(NA, t_decision[-.N]), by=trustee]

#getting rid of uncoded trustees
b = b[b$trustee > 0]
#getting rid of missed responses
#b = b[b$s_decision != 0]
#re-coding missed reponses
b$missed = 0
b$missed[b$s_decision == 0] = 1


#recoding variables
b$t_decision[b$t_decision==-1]=0
b$s_decision[b$s_decision==-1]=0
b$pt_decision[b$pt_decision==-1]=0



#factors
b$trustee = as.factor(b$trustee)
b$subject = as.factor(b$subject)
b$t_decision = as.factor(b$t_decision)
b$pt_decision = as.factor(b$pt_decision)
b$missed = as.factor(b$missed)

# add previous subject decision
b[, s_decision_lag1:=c(NA, s_decision[-.N]), by=trustee]
b[, s_decision_lag2:=c(NA, s_decision_lag1[-.N]), by=trustee]
b[, s_decision_lag3:=c(NA, s_decision_lag2[-.N]), by=trustee]
b[, s_decision_lag4:=c(NA, s_decision_lag3[-.N]), by=trustee]


#recalculating RT
b$RT_calc = b$feedback_Onset-b$decision_Onset

#mean-centering exchange
b$mc_exchange = b$exchange - 24.5

#calculating aggregate means
#RTmeansBYsubjectBYtrustee=aggregate(b[,b$RT_calc], list(b$subject,b$trustee),mean)
RTmeansBYsubject=aggregate(b[,b$RT_calc], list(b$subject),mean)
View(RTmeansBYsubject)

attach(b)
#running glmer: w/ varying intercept and slope for subject
m0 <- glmer(s_decision ~ trustee*mc_exchange*pt_decision + (1|subject), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m0)
summary(m0)
#for pt_decision and subject and varying slope effect of trustee within subject
m1 <- glmer(s_decision ~ trustee*mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
m2 <- glmer(s_decision ~ trustee*mc_exchange+trustee*pt_decision+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m1,m2)

#m3=best model for beha data#
m3.beh <- glmer(s_decision ~ trustee*pt_decision+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
m4.beh <- glmer(s_decision ~ trustee+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m3.beh, m4.beh)
m5.beh <- glmer(s_decision ~ trustee*pt_decision+mc_exchange*pt_decision + (1|subject), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m5.beh, m3.beh)

##the best model for scan data##
m3 <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m2,m3)
anova(m3)
summary(m3)

m4 <- glmer(s_decision ~ trustee+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m3,m4)

#establishing a reference grid
likelihood.m3 <- ref.grid(m3)

#marginal means of trustree, plus contrast
trustee.lsm <- lsmeans(likelihood.m3, "trustee")
trustee.sum <- summary(trustee.lsm, infer = c(TRUE,TRUE), level = 0.95, adjust = "bon")
contrast(trustee.lsm, method = "pairwise", adjust ="bon")

# plot main effect of trustee
trustee.plot1 <- plot(lsmeans(likelihood.m3, "trustee"))
trustee.plot2 <- plot(trustee.lsm, type ~ trustee, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Trustee Type")
#axis(side = "bottom",labels = c("good","bad","neutral","computer"), at=1:4)

#rg2trust <- lsmeans(m3,by = "trustee")
#pairs(rg2trust)
#plot(rg2trust, type ~ trustee, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Trustee Type")

# plot trustee*mc_exchange
trusteeXexchange.lsm <- lsmeans(likelihood.m3, "mc_exchange", by = "trustee")
trusteeXexchange.sum <- summary(trusteeXexchange.lsm, infer = c(TRUE,TRUE), level = 0.95, adjust = "bon")

rg2exch <- lsmeans(m3,"mc_exchange", by = "trustee", at = list(mc_exchange = c(-24.5,  0,  24.5)))
plot(rg2exch, type ~ trustee, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Number of Exchanges")
pairs(rg2exch)
rg2exchB <- lsmeans(m3,"mc_exchange", by = "trustee", at = list(mc_exchange = c(-24.5, 24.5)))
pairs(rg2exchB)

# plot pt_decision*mc_exchange
rg3exch <- lsmeans(m3,"pt_decision", by = "mc_exchange", at = list(mc_exchange = c(-24.5,  0,  24.5)))
pt_decisionXmc_exchange.plot <- plot(rg3exch, type ~ pt_decision, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Previous trustee decision")
pairs(rg3exch)
pt_decisionXmc_exchange.plot$lwd=2
pt_decisionXmc_exchange.plot$lwd

## test pt_decision*ps_decision
# first, test model with lags

m3lag <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + pt_decision*s_decision_lag1 + s_decision_lag2 + s_decision_lag3 + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m3lag)
summary(m3lag)

anova(m3,m3lag)

## ad: this looks currently (1/19/17) like the best-fitting model
# lag means that lagged choice is included and t means that trustee is nested within subject on the random side
m3lagt <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + pt_decision*s_decision_lag1 + (1+pt_decision|subject) + (1|subject:trustee), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m3lagt)
summary(m3lagt)

##analyzing RTs##
b$RT_calc_ln = log(b$RT_calc)
boxplot(split(b$RT_calc_ln,b$subject))

attach(b)
m0.RT <- lme(RT_calc_ln ~ trustee*mc_exchange*pt_decision*s_decision, random = (1|subject), data = b,na.action = na.omit)
m1.RT <- lme(RT_calc_ln ~ trustee*pt_decision*mc_exchange*s_decision + (1+pt_decision|subject), data = b,na.action = na.omit)
m2.RT <- lme(RT_calc_ln ~ trustee*pt_decision*mc_exchange*s_decision + (1+pt_decision|subject)+(1+s_decision|subject), data = b,na.action = na.omit)
m3.RT <- lme(RT_calc_ln ~ trustee*pt_decision*mc_exchange*s_decision + (1+s_decision|subject), data = b,na.action = na.omit)

##analyzing missed responses###
missedBYsubject=aggregate(b[,b$missed], list(b$subject),sum)
m0.missed <-glmer(missed ~ trustee*mc_exchange + (1|subject), binomial(link = "logit"), data = b,na.action = na.omit)
