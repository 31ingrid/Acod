
library(batch)

Allk<-c(rep(218,50))  
SSBA<-c(10914)#,2728,5457,10914)
seed = 111

for(i in 1:length(Allk)) 
for (j in 1:length(SSBA))
{
{
 seed <- rbatch("platinum_hisel.r",seed=seed,mig=Allk[i],ssba=SSBA[j])  # add whatever other arguments you want to pass to the AveRes.R function inside the rbatch command
}
}



