
#library 
library("R.matlab")
library("reshape2")
library("lme4")
library("ggplot2")
library("gridExtra")
library("MuMIn")
library("corrplot")
library("lmerTest")
library("lmtest")


setwd("C:/Users/feldhaus/Desktop/Dokumente/05_Diss/03_calc_study2/")


#1 load data 
#2 mixed model testing for an influence of the treatment protocols
#3 ANOVA testing for a difference relief in depending on treatment protocol and session
#4 mixed model for session 1 and session 2
#5 regression analysis testing for the influence of the different blocks (start, end)

#############################################################################################

#### 1 load data 

#load data
temp <- readMat("mean_ratings.mat")
mean_ratings <- temp$mean.ratings
temp <- readMat("Ind_start.mat")
Ind_start <- temp$inc.dec
temp <- readMat("mean_diffs.mat")
mean_diffs <- temp$mean.diffs

# ordered by frequencies (50Hz, 100Hz, 150Hz) - not by time
diffs_sess1 <- data.frame(read.csv("diffs_s1.csv"))
diffs_sess2_wo <- data.frame(read.csv("diffs_s2.csv")) #wrong order
diffs_sess2 <- data.frame(rbind(diffs_sess2_wo[which(diffs_sess2_wo$start == 2),], 
                                diffs_sess2_wo[which(diffs_sess2_wo$start == 1),]))

diffs_both_org <- rbind(cbind(diffs_sess1, SbId = seq(11:31)), cbind(diffs_sess2, SbId = seq(11:31)))
diffs_both_org$sess <- c(rep(1,21), rep(2,21))


# now it will be ordered by time instead of frequency!
diffs_both <- data.frame(matrix(nrow = 42, ncol = 7)) # SbId, 8 ratings, start indice
diffs_both[,1] <- c(seq(11,31,by = 2), seq(12,30, by = 2))
diffs_both[,6] <- c(diffs_sess1[,5], diffs_sess2[,5])
diffs_both[,5] <- diffs_both_org[,4]
diffs_both[,7] <- c(rep(1,dim(diffs_sess1)[1]), rep(2,dim(diffs_sess1)[1]))
for (i in 1:dim(diffs_both)[1]){
  if (diffs_both[i,6] == 1){
    diffs_both[i,2:4] <- diffs_both_org[i,c(3,2,1)]
  } else {
    diffs_both[i,2:4] <- diffs_both_org[i,c(1,2,3)]
  }
}
colnames(diffs_both) <- c("SbId","diff1", "diff2", "diff3", "diff4", "start", "session")


#############################################################################################

#### 2 mixed model testing for an influence of the treatment protocols

#all participants in one 
#Mixed model
mean_rats1 <- rbind(cbind(mean_ratings[,c(1,8,9)],rep(1,dim(mean_ratings)[1])),cbind(mean_ratings[,c(1,18,19)],rep(2,dim(mean_ratings)[1])))
colnames(mean_rats1) <- c("SbId", "ctrl4", "trt4", "start")
mean_rats_lf <- melt(data = data.frame(mean_rats1), id.vars = c("SbId", "start"), measure.vars = c("ctrl4", "trt4"))
mean_rats_lf2 <- cbind(mean_rats_lf, cond2 = rep(c(1,2), each = dim(mean_rats1)[1]))
mean_rats_lf3 <- cbind(mean_rats_lf2, SbId2 = paste0(mean_rats_lf2$SbId,mean_rats_lf2$start))

model.null <- lmer(value ~ variable + 
                     (variable|SbId), 
                   data = mean_rats_lf3)
model.null2 <- lmer(value ~ variable + start + 
                     (variable|SbId), 
                   data = mean_rats_lf3)
model.one <- lmer(value ~ variable*start +
                    (variable|SbId), 
                  data = mean_rats_lf3)

r.squaredGLMM(model.null2)
r.squaredGLMM(model.one)
anova(model.null, model.null2)
mc1 <- anova(model.null2, model.one, refit = F)
aov(model.null2, model.one)
summary(model.null)
m02 <- summary(model.null2)
m1 <- summary(model.one)
#lrtest(model.null, model.null2)
#lrtest(model.null2, model.one)

g1 <- ggplot(mean_rats_lf2, aes(x = cond2, y = value, colour = as.factor(start))) +
  facet_wrap(~SbId, nrow=3) +
  geom_point() +
  theme_classic() +
  geom_line(data = cbind(mean_rats_lf2, pred = predict(model.null)), 
            aes(x = cond2, y = pred, colour = as.factor(start)), size = 1) +
  theme(legend.position = "none") + ylab("Pain rating (VAS)") +
  scale_color_manual(values=c(rep(c('#a2142fff','#0072bdff'))))+
  scale_x_continuous(labels=c("1" = "Ctrl", "2" = "Trt"), breaks = c(1,2), limits = c(0,3)) 
g1b <- ggplot(mean_rats_lf2, aes(x = cond2, y = value, colour = as.factor(start))) +
  facet_wrap(~SbId, nrow=3) +
  geom_point() +
  theme_classic() +
  geom_line(data = cbind(mean_rats_lf2, pred = predict(model.null2)), 
            aes(x = cond2, y = pred, colour = as.factor(start)), size = 1) +
  theme(legend.position = "none") + ylab("Pain rating (VAS)") +
  scale_color_manual(values=c(rep(c('#a2142fff','#0072bdff'))))+
  scale_x_continuous(labels=c("1" = "Ctrl", "2" = "Trt"), breaks = c(1,2), limits = c(0,3)) 
g2 <- ggplot(mean_rats_lf2, aes(x = cond2, y = value, colour = as.factor(start))) +
  facet_wrap(~SbId, nrow=3) +
  geom_point() +
  theme_classic() +
  geom_line(data = cbind(mean_rats_lf2, pred = predict(model.one)), 
            aes(x = cond2, y = pred, colour = as.factor(start)), size = 1) +
  theme(legend.position = "none") + ylab("Pain rating (VAS)") +
  scale_color_manual(values=c(rep(c('#a2142fff','#0072bdff'))))+
  scale_x_continuous(labels=c("1" = "Ctrl", "2" = "Trt"), breaks = c(1,2), limits = c(0,3)) 


#############################################################################################

#### 3 ANOVA testing for a difference relief in depending on treatment protocol and session

#ANOVA testing for an interaction of start and session
both_lf <- melt(data = diffs_both, id.vars = c("SbId", "start", "session","diff1", "diff2", "diff3"), 
                measure.vars = c("diff4"))
colnames(both_lf)[8] <- c("relief") 
aov_m <- aov(relief ~ session * start, both_lf)
summary(aov_m)
# post-hoc t-tests
t.test(both_lf[which(both_lf$start == 1 & both_lf$session == 1),8], 
       both_lf[which(both_lf$start == 2 & both_lf$session == 1),8])
t.test(both_lf[which(both_lf$start == 1 & both_lf$session == 2),8], 
       both_lf[which(both_lf$start == 2 & both_lf$session == 2),8])


aov_m2 <- aov(diff4 ~ sess, diffs_both_org)
summary(aov_m2)

se <- function(x){sd(x, na.rm = T)/length(x)}

both_lf_m_se <- data.frame(matrix(NA, 4, 4))
colnames(both_lf_m_se) <- c("start", "session", "relief", "se")
starti <- 1
for (i in 1:2){
  for (j in 1:2){
    both_lf_m_se[starti,1] <- i
    both_lf_m_se[starti,2] <- j
    both_lf_m_se[starti,3] <- mean(both_lf[which((both_lf$start == i) & (both_lf$session == j)),8])
    both_lf_m_se[starti,4] <- se(both_lf[which((both_lf$start == i) & (both_lf$session == j)),8])
    starti <- starti + 1
  }
}
both_lf_m_se$start   <- factor(both_lf_m_se$start)
both_lf_m_se$session <- factor(both_lf_m_se$session)


ggplot(both_lf_m_se, aes(x=session, y = relief, group = start, color = start, fill = start)) + 
  geom_line(size = 1) +
  geom_pointrange(aes(ymin=relief-se, ymax=relief+se), fatten = 2, size = 1) + 
  scale_color_manual(values=c(rep(c('#0072bdff','#a2142fff'))), 
                     name = "", labels = c("strong start","weak start")) +
  theme_classic() + xlab("Session") + ylab("Pain relief")

both_lf_m_se$session_jit <- c(0.8,1.8,1.2,2.2)

ggplot(both_lf_m_se, aes(x=session_jit, y = relief, group = start, color = start, fill = start)) + 
  geom_line(size = 1) +
  geom_bar(stat="identity", width = 0.18) +
  geom_pointrange(aes(ymin=relief-se, ymax=relief+se), fatten = 2, size = 1) + 
  scale_fill_manual(values=alpha(c(rep(c('#0072bdff','#a2142fff'))),0.4), 
                    name = "", labels = c("strong start","weak start")) +
  scale_color_manual(values=c(rep(c('#0072bdff','#a2142fff'))), 
                     name = "", labels = c("strong start","weak start")) +
  theme_classic() + xlab("Session") + ylab("Treatment relief 100Hz test") + 
  scale_x_continuous(breaks = c(1,2), limits = c(0.5,2.5))


#############################################################################################

#### 4 mixed model for session 1 and session 2

#devide in session1 session2
mean_ratings_sess1 <- matrix(nrow = 21, ncol = 10) # SbId, 8 ratings, start indice
mean_ratings_sess1[,1] <- mean_ratings[,1]
mean_ratings_sess1[,10] <- Ind_start[,2]
for (i in 1:dim(mean_ratings)[1]){
  if (Ind_start[i,2] == 2){
    mean_ratings_sess1[i,2:9] <- mean_ratings[i,12:19]
  } else {
    mean_ratings_sess1[i,2:9] <- mean_ratings[i,2:9]
  }
}
colnames(mean_ratings_sess1) <- c("SbId", "ctrl1", "trt1", "ctrl2", "trt2", "ctrl3", "trt3", "ctrl4", "trt4", "start")

#Mixed model
mean_rats_sess1_lf <- melt(data = data.frame(mean_ratings_sess1), id.vars = c("SbId", "start"), measure.vars = c("ctrl4", "trt4"))
mean_rats_sess1_lf[,5] <- cbind(rep(c(1,2),each=21))
model.nulla <- lmer(value ~ variable + 
                     (1|SbId), 
                   data = mean_rats_sess1_lf)
model.nulla2 <- lmer(value ~ variable + start +
                      (1|SbId), 
                    data = mean_rats_sess1_lf)
model.onea <- lmer(value ~ variable*start +
                    (1|SbId), 
                  data = mean_rats_sess1_lf)
r_n <- r.squaredGLMM(model.nulla)
r_n2 <- r.squaredGLMM(model.nulla2)
r_o <- r.squaredGLMM(model.onea)
paste0("Rconditional between nullmodel and fullmodel", round(r_o[2],3), " - ",round(r_n2[2],3), " = ",round((r_o[2]-r_n2[2]),3))

summary(model.nulla)
m02a <- summary(model.nulla2)
m1a  <- summary(model.onea)
mc2 <- anova(model.nulla2, model.onea, refit = F)
lrtest(model.nulla, model.nulla2)
lrtest(model.nulla2, model.onea)

g3 <- ggplot(mean_rats_sess1_lf, aes(x = V5, y = value, colour = as.factor(SbId))) +
  geom_point() +
  theme_classic() +
  geom_line(data = cbind(mean_rats_sess1_lf, pred = predict(model.null)), 
            aes(x = V5, y = pred, colour = as.factor(SbId))) +
  scale_color_manual(values=c(rep(c('#0072bdff','#a2142fff'),times = 10),'#0072bdff'))+
  ylab("Pain rating (VAS)") +
  theme(legend.position = "none") +
  scale_x_continuous(labels=c("1" = "Ctrl", "2" = "Trt"), breaks = c(1,2), limits = c(0,3), name = "Condition") 

g3b <- ggplot(mean_rats_sess1_lf, aes(x = V5, y = value, colour = as.factor(SbId))) +
  geom_point() +
  theme_classic() +
  geom_line(data = cbind(mean_rats_sess1_lf, pred = predict(model.null2)), 
            aes(x = V5, y = pred, colour = as.factor(SbId))) +
  scale_color_manual(values=c(rep(c('#0072bdff','#a2142fff'),times = 10),'#0072bdff'))+
  ylab("Pain rating (VAS)") +
  theme(legend.position = "none") +
  scale_x_continuous(labels=c("1" = "Ctrl", "2" = "Trt"), breaks = c(1,2), limits = c(0,3), name = "Condition") 

g4 <- ggplot(mean_rats_sess1_lf, aes(x = V5, y = value, colour = as.factor(SbId))) +
  geom_point() +
  theme_classic() +
  geom_line(data = cbind(mean_rats_sess1_lf, pred = predict(model.one)), 
            aes(x = V5, y = pred, colour = as.factor(SbId))) +
  scale_color_manual(values=c(rep(c('#0072bdff','#a2142fff'),times = 10),'#0072bdff'))+
  ylab("Pain rating (VAS)") +
  theme(legend.position = "none") +
  scale_x_continuous(labels=c("1" = "Ctrl", "2" = "Trt"), breaks = c(1,2), limits = c(0,3), name = "Condition") 

#devide in session1 session2
mean_ratings_sess2 <- matrix(nrow = 21, ncol = 10) # SbId, 8 ratings, start indice
mean_ratings_sess2[,1] <- mean_ratings[,1]
mean_ratings_sess2[,10] <- Ind_start[,2]
for (i in 1:dim(mean_ratings)[1]){
  if (Ind_start[i,2] == 1){
    mean_ratings_sess2[i,2:9] <- mean_ratings[i,12:19]
  } else {
    mean_ratings_sess2[i,2:9] <- mean_ratings[i,2:9]
  }
}
colnames(mean_ratings_sess2) <- c("SbId", "ctrl1", "trt1", "ctrl2", "trt2", "ctrl3", "trt3", "ctrl4", "trt4", "start")

#Mixed model
mean_rats_sess2_lf <- melt(data = data.frame(mean_ratings_sess2), id.vars = c("SbId", "start"), measure.vars = c("ctrl4", "trt4"))
mean_rats_sess2_lf[,5] <- cbind(rep(c(1,2),each=21))
model.nullb <- lmer(value ~ variable + 
                      (1|SbId), 
                    data = mean_rats_sess2_lf)
model.nullb2 <- lmer(value ~ variable + start +
                      (1|SbId), 
                    data = mean_rats_sess2_lf)
model.oneb <- lmer(value ~ variable*start +
                     (1|SbId), 
                   data = mean_rats_sess2_lf)

anova(model.nullb,model.oneb)
anova(model.nullb,model.nullb2)
mc3 <- anova(model.nullb2,model.oneb, refit = F)
summary(model.nullb)
m02b <- summary(model.nullb2)
m1b  <- summary(model.oneb)
lrtest(model.nullb, model.nullb2)
lrtest(model.nullb2, model.oneb)

r_nb2 <- r.squaredGLMM(model.nullb2)
r_ob <- r.squaredGLMM(model.oneb)
paste0("Rconditional between nullmodel and fullmodel", round(r_ob[2],3), " - ",round(r_nb2[2],3), " = ",round((r_ob[2]-r_nb2[2]),3))

g5 <- ggplot(mean_rats_sess2_lf, aes(x = V5, y = value, colour = as.factor(SbId))) +
  geom_point() +
  theme_classic() +
  geom_line(data = cbind(mean_rats_sess2_lf, pred = predict(model.nullb)), 
            aes(x = V5, y = pred, colour = as.factor(SbId))) +
  scale_color_manual(values=c(rep(c('#0072bdff','#a2142fff'),times = 10),'#0072bdff'))+
  ylab("Pain rating (VAS)") +
  theme(legend.position = "none") +
  scale_x_continuous(labels=c("1" = "Ctrl", "2" = "Trt"), breaks = c(1,2), limits = c(0,3), name = "Condition") 

g5b <- ggplot(mean_rats_sess2_lf, aes(x = V5, y = value, colour = as.factor(SbId))) +
  geom_point() +
  theme_classic() +
  geom_line(data = cbind(mean_rats_sess2_lf, pred = predict(model.nullb2)), 
            aes(x = V5, y = pred, colour = as.factor(SbId))) +
  scale_color_manual(values=c(rep(c('#0072bdff','#a2142fff'),times = 10),'#0072bdff'))+
  ylab("Pain rating (VAS)") +
  theme(legend.position = "none") +
  scale_x_continuous(labels=c("1" = "Ctrl", "2" = "Trt"), breaks = c(1,2), limits = c(0,3), name = "Condition") 


g6 <- ggplot(mean_rats_sess2_lf, aes(x = V5, y = value, colour = as.factor(SbId))) +
  geom_point() +
  theme_classic() +
  geom_line(data = cbind(mean_rats_sess2_lf, pred = predict(model.oneb)), 
            aes(x = V5, y = pred, colour = as.factor(SbId))) +
  scale_color_manual(values=c(rep(c('#0072bdff','#a2142fff'),times = 10),'#0072bdff'))+
  ylab("Pain rating (VAS)") +
  theme(legend.position = "none") +
  scale_x_continuous(labels=c("1" = "Ctrl", "2" = "Trt"), breaks = c(1,2), limits = c(0,3), name = "Condition") 

grid.arrange(g3b,g4,g5b,g6, nrow = 2)

sum_models <- matrix(NA, 26, 6)
sum_models[3,2:6] <- round(m02$coefficients[2,],3)
sum_models[4,2:6] <- round(m02$coefficients[3,],3)
sum_models[6,2:6] <- round(m1$coefficients[2,],3)
sum_models[7,2:6] <- round(m1$coefficients[3,],3)
sum_models[8,2:6] <- round(m1$coefficients[4,],3)

sum_models[12,2:6] <- round(m02a$coefficients[2,],3)
sum_models[13,2:6] <- round(m02a$coefficients[3,],3)
sum_models[15,2:6] <- round(m1a$coefficients[2,],3)
sum_models[16,2:6] <- round(m1a$coefficients[3,],3)
sum_models[17,2:6] <- round(m1a$coefficients[4,],3)

sum_models[21,2:6] <- round(m02b$coefficients[2,],3)
sum_models[22,2:6] <- round(m02b$coefficients[3,],3)
sum_models[24,2:6] <- round(m1b$coefficients[2,],3)
sum_models[25,2:6] <- round(m1b$coefficients[3,],3)
sum_models[26,2:6] <- round(m1b$coefficients[4,],3)

mod_comp <- matrix(NA, 3, 5)
for (i in 1:3){
  if (i == 1){ 
    mc = mc1
  } if (i == 2){
    mc = mc2
  } if (i == 3){
    mc = mc3
  }
  mod_comp[i,2] = mc$Df
}
mod_comp[]

#############################################################################################

#### 5 regression analysis testing for the influence of the different blocks (start, end)

fit_both <- lm(diff4 ~ diff1 + diff2 + diff3, diffs_both)
summary(fit_both) # show results

gx <- ggplot(data = diffs_both) +
  geom_point(aes(x = diff1, y = diff4), color = "lightblue") +
  geom_point(aes(x = diff2, y = diff4), color = "black") +
  geom_point(aes(x = diff3, y = diff4), color = "green") +
  geom_smooth(aes(x = diff1, y = diff4), formula = y~x, color = "lightblue", method=lm, se = F) +
  geom_smooth(aes(x = diff2, y = diff4), formula = y~x, color = "black", method=lm, se = F) +
  geom_smooth(aes(x = diff3, y = diff4), formula = y~x, color = "green", method=lm, se = F) +
  xlab("Pain relief pre") + ylab("Pain relief test") +
  theme_classic()

sess1 <- diffs_both[c(1:21),]  
fit_1 <- lm(diff4 ~ diff1 + diff3, sess1)
summary(fit_1) # show results

fit_1b <- lm(diff4 ~ diff3 + diff1, sess1)
summary(fit_1b) # show results

gy

ggplot(data = sess1) +
  geom_point(aes(x = diff1, y = diff4), color = "purple", size = 2, alpha = .6) +
  #geom_point(aes(x = diff2, y = diff4), color = "black") +
  geom_point(aes(x = diff3, y = diff4), color = "green", size = 2, alpha = .6) +
  geom_smooth(aes(x = diff1, y = diff4), formula = y~x, color = "purple", method=lm, se = F) +
  #geom_smooth(aes(x = diff2, y = diff4), formula = y~x, color = "black", method=lm, se = F) +
  geom_smooth(aes(x = diff3, y = diff4), formula = y~x, color = "green", method=lm, se = F) +
  xlab("Pain relief pre") + ylab("Pain relief test") +
  theme_classic()
 
sess1_lf <- melt(data = sess1, measure.vars = c("diff1", "diff3"), 
                 id.vars = c("SbId", "diff2", "diff4", "start", "session"))
sess1_lf$col <- c(rep(1,10), rep(2,11), rep(3,10), rep(4,11))
sess1_lf1 <- sess1_lf[1:21,]

ggplot(data = sess1_lf, aes(x = value, y = diff4, color = paste(variable))) +
  geom_point(size = 2) +
  scale_color_manual(values=alpha(c(rep(c("green", "purple"))),1), 
                    name = "", labels = c("1st cond","3rd cond")) +
  geom_smooth(data = sess1, aes(x = diff1, y = diff4), formula = y~x, color = "green", method=lm, se = F) +
  geom_smooth(data = sess1, aes(x = diff3, y = diff4), formula = y~x, color = "purple", method=lm, se = F) +
  xlab("Pain relief pre") + ylab("Pain relief test") +
  theme_classic()

sess2 <- diffs_both[c(22:42),]  
fit_2 <- lm(diff4 ~ diff1 + diff2 + diff3, sess2)
summary(fit_2) # show results

gz <- ggplot(data = sess2) +
  geom_point(aes(x = diff1, y = diff4), color = "lightblue") +
  geom_point(aes(x = diff2, y = diff4), color = "black") +
  geom_point(aes(x = diff3, y = diff4), color = "green") +
  geom_smooth(aes(x = diff1, y = diff4), formula = y~x, color = "lightblue", method=lm, se = F) +
  geom_smooth(aes(x = diff2, y = diff4), formula = y~x, color = "black", method=lm, se = F) +
  geom_smooth(aes(x = diff3, y = diff4), formula = y~x, color = "green", method=lm, se = F) +
  xlab("Pain relief pre") + ylab("Pain relief test") +
  theme_classic()

grid.arrange(gx,gy,gz, nrow = 1)

cbind(sess1,sess2)

sess1sess2 = cbind(sess1[,2:5],sess2[,2:5])
colnames(sess1sess2) = c("s1diff1","s1diff2","s1diff3","s1diff4",
                         "s2diff1","s2diff2","s2diff3","s2diff4")
corrplot(cor(sess1sess2), method = "ellipse")
