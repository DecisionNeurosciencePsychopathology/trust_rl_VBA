library(lme4)
library(lsmeans)
library(data.table)
library(nlme)
library(xtable)
library(ggplot2)
# load data
library(readr)

b = read_delim("~/Google Drive/skinner/projects_analyses/Project Trust/data/compiled_betas_NullTrusteeHybridHRegret_n19_RPE_mask.csv",
                            "\t", escape_double = FALSE, trim_ws = TRUE)
#beha_behavior = read_delim("~/Box Sync/Project Trust Game/data/processed/scan_behavior.txt",
#                           "\t", escape_double = FALSE, trim_ws = TRUE)
#b = data.table(beha_behavior)
b = data.table(b)
#attach(b)
View(b)

#factors
b$Model = as.factor(b$Model)
b$subject = as.factor(b$subject)
b$striatum <- b$PEs1

attach(b)
#running glmer: w/ varying intercept and slope for subject
m1 <- lmer(striatum ~ Model + (1|subject), data = b)
summary(m1)
car::Anova(m1)

colnames(b)[3] <- "VS"
colnames(b)[4] <- "a. insula"
colnames(b)[5] <- "midbrain"
colnames(b)[6] <- "DS"
colnames(b)[7] <- "visual cortex"

View(b)



## try the region as factor
bt <-  b[ ,-c(8)]
library(reshape2)
br <- melt(bt,id = c("subject","Model"))
colnames(br)[3] <- "region"
colnames(br)[4] <- "beta"

View(br)

m2 <- lmer(beta ~ Model + region +  (1|subject), data = br)
summary(m2)
car::Anova(m2)

ls_m2 <- lsmeans(m2,"Model", by = "region")

contrast(ls_m2, method = "eff", adjust ="tukey")

setwd("~/Google Drive/skinner/projects_analyses/Project Trust/data/")
pdf("ROI model comparison.pdf", width=10, height=6)
plot(ls_m2, type ~ striatum, horiz=F,ylab = "Response", xlab = "Model", color)
dev.off()

sink("ROI_model_comparison.txt")
summary(m2)
car::Anova(m2)
ls_m2 <- lsmeans(m2,"Model", by = "region")
contrast(ls_m2, method = "eff", adjust ="tukey")
sink()

library(multcomp)
m3 <- lmer(beta ~ Model*region +  (1|subject:region) + (1|subject), data = br)
summary(m3)
car::Anova(m3)
anova(m2,m3)


ls_m3 <- lsmeans(m3,"Model", by = "region")
cld(ls_m3)
contrast(ls_m3, method = "eff", adjust ="tukey", interaction = FALSE)
plot(ls_m3, type ~ striatum, horiz=F,ylab = "Response", xlab = "Model")

m4 <- lmer(beta ~ Model + region +  (1|subject:region) + (1|subject), data = br)
sink("ROI_model_comparison_nested.txt")
summary(m4)
car::Anova(m4)
ls_m4 <- lsmeans(m4,"Model")
contrast(ls_m4, method = "eff", adjust ="tukey", interaction = TRUE)
sink()

cld(ls_m4)

pdf("ROI model comparison nested.pdf", width=10, height=6)
ls_m4 <- lsmeans(m4,"Model")
plot(ls_m4, type ~ response, horiz=F,ylab = "Response, A.U.", xlab = "Model")
dev.off()

ls_m1 <- lsmeans(m1,"Model")
plot(ls_m1, type ~ striatum, horiz=F,ylab = "Striatal response", xlab = "Model")
