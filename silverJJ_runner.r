

library(batch)

Allk<-c(rep(2750,50))
SSBA<-c(15000)
seed = 51
for (j in 1:length(SSBA))
for(i in 1:length(Allk))
{
{
 seed <- rbatch("silverJJ.r",seed=seed,mig=Allk[i],ssba=SSBA[j])  # add whatever other arguments you want to pass to the AveRes.R function inside the rbatch command
}
}  


