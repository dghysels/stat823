# load packages
library(knitr)
library(formatR)
library(stargazer)
library(xtable)
library(readxl)
library(ggplot2)
library(ggcorrplot)
library(gridExtra)
library(leaps)
library(car)
library(readr)
library(faraway)
knitr::opts_chunk$set(echo = FALSE,purl = TRUE)
options(digits = 5, width = 60, xtable.comment = FALSE)#,knitr.table.format="latex")
opts_chunk$set(tidy.opts = list(width.cutoff=60), tidy=TRUE)
out_type <- knitr::opts_knit$get("rmarkdown.pandoc.to")

dataraw <- read_xlsx("../823data/pcancer.xlsx")

data2 <- data.frame(dataraw)
data2$seminal <- factor(dataraw$seminal)
levels(data2$seminal) <- c("absence", "presence")
data2$score <- factor(dataraw$score)
levels(data2$score) <- c("6", "7", "8")

attach(data2)

library(epiDisplay)
summ.table <- summ(data2[2:9])$table
summ.table[1:8] <- c("PSA","VOL","WT","AGE","HYP","SEM","CAP","SCORE")
colnames(summ.table) <- c("Variable","Obs. Count","Mean","Median","SD","Min.","Max.")
xtable(summ.table,label="tab:tab.summary",caption = "Summary of statistics")

library(psych)
library(GGally)
#pairs.panels(data2[c(2:6,8)])

ggpairs(data2[c(2:6,8)],ggplot2::aes())

p1 <- ggplot(data2, aes(y=cancerv,x=1)) + 
  geom_boxplot() + 
  ggtitle("(a) Cancer Volume") + 
  labs(x="") + 
  geom_jitter(position=position_jitter(0.05)) + 
  theme(plot.title = element_text(color="blue", size=10, face="bold.italic"))

p2 <- ggplot(data2, aes(y=capsular,x=1)) + 
  geom_boxplot() + 
  ggtitle("(b) Capsular penetration") + 
  labs(x="") + 
  geom_jitter(position=position_jitter(0.05))

p3 <- ggplot(data2, aes(y=hyperplasia,x=1)) + 
  geom_boxplot() + 
  ggtitle("(c) Hyperplasia") + 
  labs(x="") + 
  geom_jitter(position=position_jitter(0.05))

p4 <- ggplot(data2, aes(y=weight,x=1)) + 
  geom_boxplot() + 
  ggtitle("(d) Weight") + 
  labs(x="") + 
  geom_jitter(position=position_jitter(0.05))

p5 <- ggplot(data2, aes(y=age,x=1)) + 
  geom_boxplot() + 
  ggtitle("(e) Age") + 
  labs(x="") + 
  geom_jitter(position=position_jitter(0.05))

grid.arrange(ncol = 2, p1,p2)
grid.arrange(ncol = 2, p3,p4)
grid.arrange(ncol = 2, p5)



par(mfrow = c(1,2))
tab1(data2$seminal, main = "Distribution of Seminal")
tab1(data2$score, main="Distribution of Score")

#par(mfrow=c(2,2))
p1 <- ggplot(data2, aes(x=cancerv)) + 
  geom_histogram(binwidth = 2)

p2 <- ggplot(data2, aes(x=log(cancerv))) + 
  geom_histogram(binwidth = 1)

p3 <- ggplot(data2, aes(x=capsular)) + 
  geom_histogram(binwidth = 0.5)

grid.arrange(ncol = 2, p1,p2,p3)

lm.cancerv <- lm(data = data2, psa ~ cancerv)
lm.cancerv.log <- lm(data = data2, log(psa) ~ log(cancerv))
lm.hyperplasia <- lm(data = data2, psa ~ hyperplasia)
lm.capsular <- lm(data = data2, psa ~ capsular)
lm.weight <- lm(data = data2, psa ~ weight)
lm.weight.log <- lm(data = data2, log(psa) ~ log(weight))
lm.age <- lm(data = data2, psa ~ age)
lm.score6 <- lm(data = data2, psa ~ (score==6))
lm.score7 <- lm(data = data2, psa ~ (score==7))
lm.score8 <- lm(data = data2, psa ~ (score==8))
lm.seminal <- lm(data = data2, psa ~ seminal)

summ.cancerv <- summary(lm.cancerv)
summ.hyperplasia <- summary(lm.hyperplasia)
summ.capsular <- summary(lm.capsular)
summ.weight <- summary(lm.weight)
summ.age <- summary(lm.age)
summ.score8 <- summary(lm.score8)
summ.seminal <- summary(lm.seminal)

summ.cancerv.v <- c(summ.cancerv$coefficients[c(2,4,6,8)],summ.cancerv$adj.r.squared,summ.cancerv$fstatistic[1])
summ.hyperplasia.v <- c(summ.hyperplasia$coefficients[c(2,4,6,8)],summ.hyperplasia$adj.r.squared,summ.hyperplasia$fstatistic[1])
summ.capsular.v <- c(summ.capsular$coefficients[c(2,4,6,8)],summ.capsular$adj.r.squared,summ.capsular$fstatistic[1])
summ.age.v <- c(summ.age$coefficients[c(2,4,6,8)],summ.age$adj.r.squared,summ.age$fstatistic[1])
summ.weight.v <- c(summ.weight$coefficients[c(2,4,6,8)],summ.weight$adj.r.squared,summ.weight$fstatistic[1])
summ.score8.v <- c(summ.score8$coefficients[c(2,4,6,8)],summ.score8$adj.r.squared,summ.score8$fstatistic[1])
summ.seminal.v <- c(summ.seminal$coefficients[c(2,4,6,8)],summ.seminal$adj.r.squared,summ.seminal$fstatistic[1])


summarydf <- data.frame("VOL"=summ.cancerv.v,
                        "HYP"=summ.hyperplasia.v,
                        "CAP"=summ.capsular.v,
                        "AGE"=summ.age.v,
                        "WT"=summ.weight.v,
                        "SCORE"=summ.score8.v,
                        "SEM"=summ.seminal.v)
rownames(summarydf) <- c("beta1","std Err","t-value","p-value","adj. R2","f-statistic")


for(i in 1:ncol(summarydf)){
  summarydf[[i]] <- sapply(summarydf[,i],function(x){if(x>0.0001){
                                                        return(format(x,digits=4))
                                                      }else{
                                                        return(format(x,digits=2))
                                                      }
                                                    })
}

knitr::kable(summarydf,label = "tab.stat.summary")
                        


par(mfrow=c(1,2))
plot(y=data2$psa,x=data2$capsular, 
     xlab="PSA", 
     ylab="CAP", 
     main="CAP Fit")
abline(lm.capsular)

plot(fitted(lm.capsular),residuals(lm.capsular), 
     xlab="residuals",
     ylab="fitted",
     main = "CAP fitted \n vs residuals")
abline(h = 0)

par(mfrow=c(1,2))
plot(y=data2$psa,x=data2$cancerv, 
     xlab="PSA", 
     ylab="VOL", 
     main="VOL Fit")
abline(lm.cancerv)

plot(fitted(lm.cancerv),residuals(lm.cancerv), 
     xlab="residuals",
     ylab="fitted",
     main = "VOL fitted vs residuals")
abline(h = 0)


par(mfrow=c(1,2))
plot(y=data2$psa,x=data2$weight, 
     xlab="PSA", 
     ylab="WT", 
     main="WT Regression Fit")
abline(lm.weight)

plot(fitted(lm.weight),residuals(lm.weight), 
     xlab="residuals",
     ylab="fitted",
     main = "WT fitted vs residuals")
abline(h = 0)


par(mfrow=c(1,2))
plot(y=data2$psa,x=data2$age, 
     xlab="PSA", 
     ylab="AGE", 
     main="AGE Regression Fit")
abline(lm.age)

plot(fitted(lm.age),residuals(lm.age), 
     xlab="residuals",
     ylab="fitted",
     main = "AGE fitted vs residuals")
abline(h = 0)


par(mfrow=c(1,2))
plot(y=data2$psa,x=data2$age, 
     xlab="PSA", 
     ylab="HYP", 
     main="HYP Regression Fit")
abline(lm.age)

plot(fitted(lm.age),residuals(lm.age), 
     xlab="residuals",
     ylab="fitted",
     main = "HYP fitted vs residuals")
abline(h = 0)


ncvp <-ncvTest(lm.cancerv)$p
swp <- shapiro.test(residuals(lm.cancerv))$p.value


par(mfrow=c(1,2))
plot(y=log(data2$psa),x=log(data2$cancerv), 
     xlab="log(PSA)", 
     ylab="log(Cancer Volume)", 
     main="Cancer Volume \n Regression Fit")

abline(lm.cancerv.log)

plot(fitted(lm.cancerv.log),residuals(lm.cancerv.log), 
     xlab="residuals",
     ylab="fitted",
     main = "Log(Cancer Vol.)\n fitted vs residuals")
abline(h = 0)

ncvp <-ncvTest(lm.cancerv.log)$p
swp <- shapiro.test(residuals(lm.cancerv.log))$p.value
t<- data.frame(a=c(swp),b=c(ncvp))
knitr::kable(t,col.names =c("Shapiro-wilks p-value","ncv test p-value"),label="tab:tab4",caption = "Shapiro-wilks and ncv test")

b <- regsubsets(log(psa) ~ log(cancerv)+seminal+capsular+score, data = data2)
regsub <- summary(b)
kable(with(regsub, cbind(which, rss, adjr2, cp, bic)), digits = 4, label = "tab.covariate.sel")

lm.multi <- lm(data=data2, log(psa) ~ log(cancerv) + score + seminal)
lm.multi.summary <- summary(lm.multi)


gi <- influence(lm.multi)
par(mfrow=c(1,2))

halfnorm(cooks.distance(lm.multi), cex = 1.5, las=1)
influencePlot(lm.multi, id = list(n = 3))
