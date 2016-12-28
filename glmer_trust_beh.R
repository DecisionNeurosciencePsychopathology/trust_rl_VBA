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

#  this is our favoritissimo model
m3 <- glmer(s_decision ~ exchange.mc*pt_decision + pt_decision*trustee + (1+pt_decision|subject), binomial(link = "logit"), data = bt,na.action = na.omit, contrasts = list(trustee = "contr.sum"))
summary(m3)

#m4 <- glmer(s_decision ~ exchange.mc*pt_decision + pt_decision*trustee + (1+pt_decision*exchange.mc|subject), binomial(link = "logit"), data = bt,na.action = na.omit, contrasts = list(trustee = "contr.sum"))
#summary(m4)



anova(m2,m3)

# plot pt_decision*exchange
rg2exch <- lsmeans(m2,"pt_decision", by = "exchange.mc", at = list(exchange.mc = c(-20,  0,  20)))
plot(rg2exch, type ~ pt_decision, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Previous trustee decision")

# plot pt_decision*trustee
rg2trust <- lsmeans(m2,"pt_decision", by = "trustee")
plot(rg2trust, type ~ pt_decision, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Previous trustee decision")
