
library(batch)

Allk<-c(rep(55,1))  
SSBA<-c(5457)#,2927,10914)
seed = 1

for(i in 1:length(Allk)) 
for (j in 1:length(SSBA))
{
{
 seed <- rbatch("silverFF.r",seed=seed,mig=Allk[i],ssba=SSBA[j])  # add whatever other arguments you want to pass to the AveRes.R function inside the rbatch command
}
}



