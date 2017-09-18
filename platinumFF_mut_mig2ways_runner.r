
library(batch)

Allk<-c(rep(5,105))  #this vector is the number of migrants
SSBA<-c(5457)#15000 for outer coast, 5457 for inner fjord
seed = 1

for(i in 1:length(Allk)) 
for (j in 1:length(SSBA))
{
{
 seed <- rbatch("platinumFF_mut_mig2ways.r",seed=seed,mig=Allk[i],ssba=SSBA[j])  # add whatever other arguments you want to pass to the AveRes.R function inside the rbatch command
}
}



