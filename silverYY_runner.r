


library(batch)

Allk<-c(rep(2500,10))
seed = 416

for(i in 1:length(Allk))
{
 seed <- rbatch("silverYY.r",seed=seed,mig=Allk[i])  # add whatever other arguments you want to pass to the AveRes.R function inside the rbatch command
}
