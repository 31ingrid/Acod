
library(batch)

Allk<-c(rep(2300,110))
SSBA<-c(15000)
seed = 1
for (j in 1:length(SSBA))
for(i in 1:length(Allk))
{
{
 seed <- rbatch("silverEE.r",seed=seed,mig=Allk[i],ssba=SSBA[j])  # add whatever other arguments you want to pass to the AveRes.R function inside the rbatch command
}
} 


