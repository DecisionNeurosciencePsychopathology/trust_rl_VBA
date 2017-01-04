library(lme4)
library("lsmeans")
#install.packages("data.table")
library("data.table")

# load data
library(readr)
beha_behavior <- read_delim("~/Box Sync/Project Trust Game/data/processed/beha_behavior.txt",
"\t", escape_double = FALSE, trim_ws = TRUE)
View(beha_behavior)

b = beha_behavior
bt <- data.table(b)
#bt[, "pt_decision":=c(NA, t_decision[-.N]), by=trustee]
# recode decisions into 0 and 1
#bt$t_decision[t_decision==-1]=0
#bt$s_decision[s_decision==-1]=0
#bt$pt_decision[b$pt_decision==-1]=0
bt$trustee = as.factor(bt$trustee)

bt$subject = as.factor(bt$subject)
bt$pt_decision = as.factor(bt$pt_decision)
bt$t_decision = as.factor(bt$t_decision)

# code subject's previous decision
bt[, "ps_decision":=c(NA, s_decision[-.N]), by=trustee]

# remove RT outliers
na.rm = TRUE
  qnt <- quantile(bt$decision_RT, probs=c(.25, .75), na.rm = na.rm)
  H <- 1.5 * IQR(bt$decision_RT, na.rm = na.rm)
  bt$decision_RT_resc <- bt$decision_RT
  bt$decision_RT_resc[bt$decision_RT < (qnt[1] - H)] <- NA
  bt$decision_RT_resc[bt$decision_RT > (qnt[2] + H)] <- NA
  bt$logRT <- log(bt$decision_RT)



bt[, "decision_RT_resclag":=c(NA, decision_RT_resc[-.N]), by=trustee]
bt[, "decision_RT_resclag2":=c(NA, decision_RT_resclag[-.N]), by=trustee]
bt[, "decision_RT_resclag3":=c(NA, decision_RT_resclag2[-.N]), by=trustee]
bt[, "decision_RT_resclag4":=c(NA, decision_RT_resclag3[-.N]), by=trustee]
bt[, "decision_RT_resclag5":=c(NA, decision_RT_resclag4[-.N]), by=trustee]

# mean-center
bt$exchange.mc <- bt$exchange - mean(bt$exchange)

#attach(b)
#bt = as_tibble((b))

# run LME

m1 <- glmer(s_decision ~ exchange.mc*pt_decision + (1|subject), binomial(link = "logit"), data = bt,na.action = na.omit)
anova(m1)
rg1 <- lsmeans(m1,"pt_decision", by = "exchange.mc", at = list(exchange.mc = c(-20,  0,  20)))
plot(rg1, type ~ pt_decision, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Previous trustee decision")

m2 <- glmer(s_decision ~ exchange.mc*pt_decision + pt_decision*trustee + (1|subject), binomial(link = "logit"), data = bt,na.action = na.omit, contrasts = list(trustee = "contr.sum"))
summary(m2)

m3 <- glmer(s_decision ~ exchange.mc*pt_decision + pt_decision*trustee + (1+pt_decision|subject), binomial(link = "logit"), data = bt,na.action = na.omit, contrasts = list(trustee = "contr.sum"))
summary(m3)

m4 <- glmer(s_decision ~ exchange.mc*pt_decision + pt_decision*trustee + (1|subject) + (1|trustee:subject), binomial(link = "logit"), data = bt,na.action = na.omit, contrasts = list(trustee = "contr.sum"))
summary(m4)

#  this is our favoritissimo model
m5 <- glmer(s_decision ~ exchange.mc*pt_decision + pt_decision*trustee + (1 + pt_decision|subject) + (1|trustee:subject), binomial(link = "logit"), data = bt,na.action = na.omit, contrasts = list(trustee = "contr.sum"))
summary(m5)

m6 <- glmer(s_decision ~ exchange.mc*pt_decision + pt_decision*trustee + ps_decision + (1 + pt_decision + ps_decision|subject) + (1|trustee:subject), binomial(link = "logit"), data = bt,na.action = na.omit, contrasts = list(trustee = "contr.sum"))
summary(m6)

anova(m5,m6)


m7 <- glmer(s_decision ~ exchange.mc*pt_decision + pt_decision*trustee + ps_decision + (1 + pt_decision*exchange.mc|subject) + (1|trustee:subject), binomial(link = "logit"), data = bt,na.action = na.omit, contrasts = list(trustee = "contr.sum"))
summary(m7)
anova(m7)

# plot pt_decision*exchange
rg4exch <- lsmeans(m4,"pt_decision", by = "exchange.mc", at = list(exchange.mc = c(-20,  0,  20)))
plot(rg4exch, type ~ pt_decision, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Previous trustee decision")

rg6exch <- lsmeans(m6,"pt_decision", by = "exchange.mc", at = list(exchange.mc = c(-20,  0,  20)))
plot(rg6exch, type ~ pt_decision, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Previous trustee decision")


# plot pt_decision*trustee
rg4trust <- lsmeans(m4,"pt_decision", by = "trustee")
plot(rg4trust, type ~ pt_decision, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Previous trustee decision")

rg6trust <- lsmeans(m6,"pt_decision", by = "trustee")
plot(rg6trust, type ~ pt_decision, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Previous trustee decision")

rg6trustmain <- lsmeans(m6,"trustee")
plot(rg6trustmain, type ~ pt_decision, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Trustee")

# RT analysis
rt_m1 = lmer(decision_RT_resc ~ exchange.mc*pt_decision + pt_decision*trustee + decision_RT_resclag + (1 + pt_decision|subject) + (1|trustee:subject),  data = bt,na.action = na.omit, contrasts = list(trustee = "contr.sum"))
summary(rt_m1)
anova(rt_m1)


rt_m2 = lme(decision_RT_resc ~ exchange.mc + trustee + pt_decision +  decision_RT_resclag, random =~  1 + exchange.mc + decision_RT_resclag |subject,  data = bt,na.action = na.omit, contrasts = list(trustee = "contr.sum"))
summary(rt_m2)
anova(rt_m2)

rt_m3 = lme(decision_RT_resc ~ trialnum + trustee*pt_decision*s_decision + decision_RT_resclag, random =~  1 + trialnum + decision_RT_resclag |subject,  data = bt,na.action = na.omit, contrasts = list(trustee = "contr.sum"))
summary(rt_m3)
anova(rt_m3)



rg_rt_m2 <- lsmeans(rt_m2,"trustee")
plot(rg_rt_m2, type ~ trustee, horiz=F, ylab = "RT (predicted)", xlab = "Trustee")

