# Acod
This repository contains code that will run simulations in the manuscript entitled
"Inferring genetic connectivity in real populations, exemplified by coastal and oceanic Atlantic cod" by Spies et al.
Each simulation was replicated 100 times with different random number seeds. 

Data in Table S1 in the Appendix was run using the following code.

The upper section of S1 "Migration from the N. Sea to the inner fjord" used the following code. 
"silverFF.r" was used in simulations for Base case, half pop. size, and double pop. size implemented by the runner file 
"silverFF_runner.r". Change the SSBA value in the runner file to 5457 for base case inner fjord sizes, to 2927 for half population 
size of the inner fjord, and to 10914 for the double population size simulations. 

"silverLL.r" is used for lower selectivity and implemented by the runner file "silverLL_runner.r".
"platinum_hisel.r" is used for higher selectivity and implemented by the runner file "platinum_hisel_runner.r". 
Change SSBA in runner file to 5457.
"silverOO.r" is used for lower fishing mortality and implemented by runner file of the same name.
"silver SS.r" is used for higher maturity and implemented by runner file of the same name.

The section of S1 entitled "Migration from teh N. Sea to the outer coast"
"silverAA.r" and the runner file of the same name is used for base case runs (change SSBA to 15000 for base case, 7400 for half 
outer population size, and 30000 for double outer population size).
"silver BBB.r" is for runs with lower selectivity and implemented by runner file of the same name.
"platinum_hisel.r" is for runs with higher selectivity and implemented by runner file of the same name. Be sure to change SSBA 
in the runner file to 15000.
"silverEE.r" is used for lower fishing mortality and implemented by runner file of the same name.
"silverJJ.r" is used for higher maturity at age and implemented by runner file of the same name.
"silverXX.r" is used for 2xNS size and implemented by runner file of the same name.
"silverYY.r" is used for 5xNS size and implemented by runner file of the same name.

Data in Table 2 and Fig. 5 come from 
"platinumFF_lnorm.r" and the runner file "platinumFF_lnorm_runner.r". These explore stochastic migration. Change SSBA to 5,457 for inner fjord runs and 15,000
for outer coast. The vector Allk is the number of migrants.

Data in Fig. 6 was generated as follows:
"platinumFF_mut_mig2ways.r" and the runner file of the same name are the files that generated boxplot B. 
Edits to the section in lines 1616-1711 were used to change from 2-way to 1-way migration to generate boxplot A. 
The same file was edited to remove lines 1588-1612 which removed migration for boxplot E.

"silverFF_noagest.r" and the runner file of the same name were used to generate boxplots C, D, and F.


