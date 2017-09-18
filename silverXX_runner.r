


library(batch)

Allk<-c(rep(1200,2))
seed = 1

for(i in 1:length(Allk))
{
 seed <- rbatch("silverXX.r",seed=seed,mig=Allk[i])  # add whatever other arguments you want to pass to the AveRes.R function inside the rbatch command
}
