library(lme4)
library(lsmeans)
library(data.table)
library(nlme)

# load data
library(readr)

scan_behavior = read_delim("~/Box Sync/Project Trust Game/data/processed/scan_behavior.txt",
                           "\t", escape_double = FALSE, trim_ws = TRUE)
beha_behavior <- read_delim("~/Box Sync/Project Trust Game/data/processed/beha_behavior.txt",
                            "\t", escape_double = FALSE, trim_ws = TRUE)
b = data.table(scan_behavior)
#b = data.table(beha_behavior)
#attach(b)

View(b)

#lag variable: previous trustee decision
b[, pt_decision:=c(NA, t_decision[-.N]), by=trustee]


#the likelihood of the next subject's decision to be share vs. keep 
#within subject's current KEEP decision (Trustee shared vs. kept)
#creating a lead variable for next subject's decision
b[, ns_decision:=shift(.SD, 1, NA, "lead"), .SDcols = "s_decision", by=subject]

#getting rid of uncoded trustees
b = b[b$trustee > 0]
#getting rid of missed responses
b = b[b$s_decision != 0]

#recoding variables
b$t_decision[b$t_decision==-1]=0
b$s_decision[b$s_decision==-1]=0
b$pt_decision[b$pt_decision==-1]=0
b$ns_decision[b$ns_decision==-1]=0

#factors
b$trustee = as.factor(b$trustee)
b$subject = as.factor(b$subject)
b$t_decision = as.factor(b$t_decision)
b$pt_decision = as.factor(b$pt_decision)
b$ns_decision = as.factor(b$ns_decision)

#recalculating RT
b$RT_calc = b$feedback_Onset-b$decision_Onset

#mean-centering exchange
b$mc_exchange = b$exchange - 24.5

#just subject decisions = keep
c = b[b$s_decision == 0]
#just trustee decisions = keep
#d = b[b$t_decision == 0]

mc0 <- glmer(ns_decision ~ t_decision + (1|subject), binomial(link = "logit"), data = c,na.action = na.omit)
mc1 <- glmer(ns_decision ~ trustee*t_decision + (1|subject), binomial(link = "logit"), data = c,na.action = na.omit)
mc2 <- glmer(ns_decision ~ t_decision*mc_exchange + (1|subject), binomial(link = "logit"), data = c,na.action = na.omit)
mc3 <- glmer(ns_decision ~ t_decision*(mc_exchange + trustee)+ (1|subject), binomial(link = "logit"), data = c,na.action = na.omit)
mc4 <- glmer(ns_decision ~ t_decision*(mc_exchange + trustee)+mc_exchange*trustee+ (1|subject), binomial(link = "logit"), data = c,na.action = na.omit)
