
library(readxl)
pcancerdf <- readxl::read_xlsx("../823data/pcancer.xlsx")

pcancerdf

pcancerdf.r <- pcancerdf[pcancerdf$weight < 400,]
pcancerdf.r <- pcancerdf.r[pcancerdf.r$psa < 100,]

plot(pcancerdf.r$psa ~ pcancerdf.r$weight)

m.p_w <- lm(log(pcancerdf.r$psa) ~ log(pcancerdf.r$weight))
summary(m.p_w)

plot(log(pcancerdf.r$psa) ~ log(pcancerdf.r$weight))
abline(m.p_w)

plot(residuals(m.p_w))

