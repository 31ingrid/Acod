

library(batch)

Allk<-c(rep(48,50),rep(49,50))
SSBA<-c(5457)
seed = 200
for (j in 1:length(SSBA))
for(i in 1:length(Allk))
{
{
 seed <- rbatch("silverSS.r",seed=seed,mig=Allk[i],ssba=SSBA[j])  # add whatever other arguments you want to pass to the AveRes.R function inside the rbatch command
}
}   




