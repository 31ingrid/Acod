
library(batch)

Allk<-c(rep(2925,155))  
SSBA<-c(15000)#,2728,5457,10914)
seed = 1

for(i in 1:length(Allk)) 
for (j in 1:length(SSBA))
{
{
 seed <- rbatch("platinumFF_lnorm.r",seed=seed,mig=Allk[i],ssba=SSBA[j])  # add whatever other arguments you want to pass to the AveRes.R function inside the rbatch command
}
}



