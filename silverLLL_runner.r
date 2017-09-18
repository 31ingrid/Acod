

library(batch)

Allk<-c(rep(73,110))
SSBA<-c(5457)
seed = 1
for (j in 1:length(SSBA))
for(i in 1:length(Allk))
{
{
 seed <- rbatch("silverLLL.r",seed=seed,mig=Allk[i],ssba=SSBA[j])  # add whatever other arguments you want to pass to the AveRes.R function inside the rbatch command
}
}    


