

library(batch)

Allk<-c(rep(1150,100),rep(2500,100))
SSBA<-c(15000)
seed = 200
for (j in 1:length(SSBA))
for(i in 1:length(Allk))
{
{
 seed <- rbatch("silverBBB.r",seed=seed,mig=Allk[i],ssba=SSBA[j])  # add whatever other arguments you want to pass to the AveRes.R function inside the rbatch command
}
}    


