

#high sel


#now make flexible with number of microsatellites
#set.seed(4);#worked
#set.seed(5);#worked
#set.seed(6);
library(batch)
parseCommandArgs()

library(Rcpp)
library(inline)
require(RcppArmadillo)
run8pops = '

using namespace Rcpp;
using namespace sugar;
RNGScope scope;

// Variables (and constants) section

NumericVector inputs(INPUTS); 

int SSBinitA; SSBinitA=inputs[0];
double pop1_prop; pop1_prop=inputs[1];
int Nyrs; Nyrs=inputs[2];
int nsamp; nsamp=inputs[3];
int ngen; ngen=inputs[4];
double fishmort; fishmort=inputs[6];
int samptype; samptype=inputs[7];
int nsims; nsims=inputs[8];
int N; N=inputs[9];
int Yrs; Yrs=inputs[10];
int mmt; mmt=inputs[11];
double dpow; dpow=inputs[12];
double rmate; rmate=inputs[13];
int ages; ages=inputs[14];
int popdyn;popdyn=inputs[15];
int nSNP; nSNP=inputs[16];
int nmig; nmig=inputs[17];
double nmig_prop;
int SSBinitB; SSBinitB=inputs[18];
int nmig1; nmig1=nmig;

NumericVector fishsel(ages); //this comes from the average selectivity at age 1963-2011
NumericVector fishselB(ages);

//low  fishsel
//fishsel[0]=0.201 //actually am using zero
fishselB[0]=0;fishselB[1]=0.656;fishselB[2]=0.891;fishselB[3]=0.783;fishselB[4]=0.742;fishselB[5]=0.743;fishselB[6]=0.743;

//high fishsel this is 65% CI 
//fishselB[0]=0.228;
//fishselB[0]=0;fishselB[1]=0.877;fishselB[2]=1;fishselB[3]=0.930;fishselB[4]=0.876;fishselB[5]=0.863;fishselB[6]=0.863;

//mean fishsel
//fishsel[0]=0.220;
fishsel[0]=0;fishsel[1]=0.818;fishsel[2]=0.974;fishsel[3]=0.895;fishsel[4]=0.847;fishsel[5]=0.848;fishsel[6]=0.848;




int npops; if (pop1_prop<1){npops=2;}else{npops=1;}
int n_l;  //proportion of fish in population 1
if (pop1_prop<0.5){n_l=1;}
if (pop1_prop==0.5){n_l=5;}
NumericVector N_til_init_a(ages); //numbers at age with 1 recruit
NumericVector N_til_init_a40(ages); //numbers at age with 1 recruit and F40%
NumericVector S_l0(2); //SSB in pop1 and 2
NumericVector R_l0A(Nyrs); //recruits in popA for each year
NumericVector R_l0B(Nyrs); //recruits in popB for each year
NumericVector R_l0A_hat(Nyrs); //recruitment with epsilon error term (same as for estimating numbers)
NumericVector R_l0B_hat(Nyrs); //
NumericVector Q_a(ages); //maturity index at age popA (Skagerrak)
NumericVector Q_b(ages); //maturity index at age popB (North Sea)

//North Sea maturity ICES
Q_b[0]=0.01;Q_b[1]=0.05;Q_b[2]=0.23;Q_b[3]=0.62;Q_b[4]=0.86;Q_b[5]=1;Q_b[6]=1;

//Skagerrak cod maturity (Knutsen 2011)
Q_a[0]= 0.03226423;Q_a[1]= 0.14982172;Q_a[2]= 0.48225745;Q_a[3]= 0.83117813;Q_a[4]= 0.96299517;Q_a[5]= 0.99278248;
Q_a[6]= 0.99862647;

double psi=1.75e-5;  //used in weight at age (von bertallanffy)
double theta=2.8571;   //used in weight at age
double Linf = 197; //wt at age
double K=0.1030;  //wt at age
NumericVector Length(ages); //length at age
NumericVector Weight_a(ages); //wt at age 
NumericVector Weight_b(ages); //wt at age
double meanFullF=0.8962449;
double varFullF=0.0232579;

NumericVector age(ages);
age=SeqLen(ages)-1;

int temp;  
double temp2;
NumericVector N_init_a(ages);
NumericVector N_init_a1(ages);  
NumericVector N_init_a2(ages);
NumericVector N_init_a1_FIRST(ages);  
NumericVector N_init_a2_FIRST(ages);  
//N_init declared below is sum of pop1_size and pop2_size (integer)
List genotypes(2);
int pop1_size;
int pop2_size;
List outmat(250);
double n_sampA;
double n_sampB;
NumericMatrix out1(nsims,Nyrs);
NumericMatrix out2(nsims,Nyrs); 
NumericMatrix out3(nsims,Nyrs);
NumericMatrix out4(nsims,Nyrs);
NumericMatrix out5(nsims,Nyrs);
NumericMatrix out6(nsims,Nyrs);
NumericMatrix out7(nsims,Nyrs);
NumericMatrix out8(nsims,Nyrs);
NumericMatrix out9(nsims,Nyrs);
NumericMatrix out10(nsims,Nyrs);
NumericMatrix out11(nsims,Nyrs);
NumericMatrix out12(nsims,Nyrs);
NumericMatrix out13(nsims,Nyrs);
NumericMatrix out14(nsims,Nyrs);
NumericMatrix out15(nsims,Nyrs);
NumericMatrix out16(nsims,Nyrs);

NumericVector out17(nsims);//1 is true Fst, 2 is initial SBPRA, 3 is initial SBPRB
NumericVector out_fst60(nsims);
NumericMatrix numsmat(nsims,ages);
NumericMatrix fst_submat(nsims,9);
NumericMatrix sig_submat(nsims,9);
NumericMatrix fst_submat60(nsims,9);
NumericMatrix sig_submat60(nsims,9);
NumericMatrix fst_submat9(nsims,9);
NumericMatrix sig_submat9(nsims,9);
NumericMatrix fst_submat99(nsims,9);
NumericMatrix sig_submat99(nsims,9);
NumericMatrix out21(nsims,Nyrs);//subsample SSB/initial SSB popA
NumericMatrix out22(nsims,Nyrs);//subsample SSB/initial SSB popB
NumericMatrix out23(nsims,Nyrs);//est SSB/SB40 popA
NumericMatrix out24(nsims,Nyrs);//est SSB/SB40 popB
NumericMatrix out25(nsims,ngen+Nyrs);//fst_vec
NumericVector genvec(nsims);


List pop1(nSNP+1); //here is where I split the data into pop1 and pop2
List pop2(nSNP+1);
List popA(Nyrs+ages);
List popB(Nyrs+ages);

int pick1;
int pick2;
double prob;
double dir;
int mover1;
int mover2;

NumericVector mean_arichA(Nyrs);
NumericVector mean_arichB(Nyrs);
double SPR_init_a;
double SPR_init_b;
double SPR_init_F40_a;
double SPR_init_F40_b;
double S_lyINIT_a;
double S_lyINIT_b;

NumericVector TACA(Nyrs);
NumericVector TACB(Nyrs);
NumericVector TAC(Nyrs);
NumericVector S_hat_syjA(Nyrs);
NumericVector S_hat_syjB(Nyrs);
NumericVector S_hat_syj(Nyrs);
NumericVector SSB_hatA(Nyrs);
NumericVector SSB_hatB(Nyrs);
NumericVector SSB_hat(Nyrs);
NumericVector SPRAtrue(Nyrs);
NumericVector SPRBtrue(Nyrs);
NumericVector FishAvec(Nyrs); //vector of optimal fishing mortality for popA (separated management)
NumericVector FishBvec(Nyrs);//vector of optimal fishing mortality for popB (separated management)

double CVa=0.6085044;
double CVs=0.779; 
double M=0.34;
double h=0.9;
double CVr=0.75;

NumericVector M2(ages);

NumericVector Surv(ages);
M2[0]=0.9;//M2[0]=1.038;
M2[1]=0.698;M2[2]=0.490;M2[3]=0.233;M2[4]=0.2;M2[5]=0.2;M2[6]=0.2;//natural mortality from ICES report
//Surv[0]=1;Surv[1]=0.43403;Surv[2]=0.18838;Surv[3]=0.08177;Surv[4]=0.03549;Surv[5]=0.01540;Surv[6]=0.00669;
Surv[0]=0.43403;Surv[1]=0.18838;Surv[2]=0.08177;Surv[3]=0.03549;Surv[4]=0.01540;Surv[5]=0.00669;Surv[6]=0.002;


NumericVector epsilonA(Nyrs);
NumericVector epsilonB(Nyrs);
NumericVector epsilon(Nyrs);
epsilonA[0]=0;
epsilonB[0]=0;
epsilon[0]=0;
int genrule;
int genrule2;
NumericVector N_hat_aA(ages);
NumericVector N_hat_aB(ages);
NumericVector N_hat_a(ages);
double sigma_ja2A=log(((1/pop1_prop)*pow(CVa,2))+1);
double sigma_ja2B=log(((1-(1/pop1_prop))*pow(CVa,2))+1);
double sigma_ja2=log(pow(CVa,2)+1);
double sigma_js2A=log(((1/pop1_prop)*pow(CVs,2))+1);
double sigma_js2B=log(((1-(1/pop1_prop))*pow(CVs,2))+1);
double sigma_js2=log(pow(CVs,2)+1);

double sigma_jr2=log(pow(CVr,2)+1);//may want to check this
double psiA;
double psiB;
double psi_star;
double rho=0.9;
double Fish=0; //set here so initially it is not a problem

double S_lyA;//true SSB
double S_lyB;
double eta_jyA;
double eta_jyB;
double eta_jy;
NumericVector temp8;
NumericVector temp9;
NumericVector temp10;
NumericVector temp11;
NumericVector temp12;
NumericVector temp13;
NumericVector temp14;
NumericVector temp15;
NumericVector temp16;
NumericVector temp17;
int temp18;
int temp19;
NumericVector Wt;
double FishA;
double FishB;
double sepruleA;
double sepruleB;
double seprule;
NumericVector Effort_tot;
NumericVector d(10); //10 is right because it is for the 10 areas
NumericVector Effort(10);
NumericVector Effort_norm;
NumericVector vec;
NumericVector spawnersA(ages);
NumericVector spawnersB(ages);
NumericVector spawnersA_rounded(ages);
NumericVector spawnersB_rounded(ages);
NumericVector spawnersA_rounded_star(ages);
NumericVector spawnersB_rounded_star(ages);
double pick_spawner_mal;
double pick_spawner_fem;
NumericVector which_spawnersA;
int some_integerA;
int some_integerB;
List tmpListA;
List tmpListB;
int pickM;
int pickF;
double pick_moverA;
double pick_moverB;

int mov_pick1;
int mov_pick2;
//NumericMatrix moverA2Brow1;
NumericMatrix RecA;
NumericMatrix RecB;
NumericMatrix recA;
NumericMatrix recB;
double rand_doub;
int samp_mom;
int samp_dad;
int rand_allele;
List tmpList_recsA(2);
List tmpList_recsB(2);

NumericVector NeAvec(Nyrs);  //Keep Ne from each year popA here
NumericVector NeBvec(Nyrs);  //Keep Ne from each year popB here


//NumericVector FSH;
//IntegerVector FSH_int;
List tmpdoneA(2);
List tmpdoneB(2);
NumericMatrix tmpageA_fem;
NumericMatrix tmpageA_mal;
NumericMatrix tmpageB_fem;
NumericMatrix tmpageB_mal;
NumericVector all_richA;
NumericVector all_richB;
NumericVector arichA(Nyrs);
NumericVector arichB(Nyrs);
NumericMatrix allA(Nyrs+1,ages);
NumericMatrix allB(Nyrs+1,ages);
NumericMatrix allspawnA(Nyrs,ages);
NumericMatrix allspawnB(Nyrs,ages);
NumericVector N_aA(ages);
NumericVector N_aB(ages);
List spawn_malA(Nyrs);
List spawn_femA(Nyrs);
List spawn_malB(Nyrs);
List spawn_femB(Nyrs);
NumericMatrix tmpMatrixA_1;
NumericMatrix tmpMatrixA_2;
NumericMatrix tmpMatrixB_1;
NumericMatrix tmpMatrixB_2;
IntegerVector NspawnersA(Nyrs);
IntegerVector NspawnersB(Nyrs);
NumericMatrix spawn_malA_thisyr; //twice as many genotypes as spawners
NumericMatrix spawn_femA_thisyr;
NumericMatrix spawn_malB_thisyr;
NumericMatrix spawn_femB_thisyr;
NumericMatrix tmpMatrixA_1plus;
NumericMatrix tmpMatrixA_2plus;
NumericMatrix tmpMatrixB_1plus;
NumericMatrix tmpMatrixB_2plus;
NumericMatrix tmpMatrixA_1details;
NumericMatrix tmpMatrixA_2details;
NumericMatrix tmpMatrixB_1details;
NumericMatrix tmpMatrixB_2details;
List tmpList_plusgroupA(2);
List tmpList_plusgroupB(2);
List tmpList_detailsA(2);
List tmpList_detailsB(2);
IntegerVector mal_spawnersA; //after random mating  **also make sure you get the right row
IntegerVector fem_spawnersA; //after random mating
IntegerVector mal_spawnersB; //after random mating  **also make sure you get the right row
IntegerVector fem_spawnersB; //after random mating
IntegerVector spawner_listM; //list of ones to choose from
IntegerVector spawner_listF;
NumericVector N_init_a1_star(ages);
NumericVector N_init_a2_star(ages);
NumericVector NtotalA(Nyrs+ages-1);
NumericVector NtotalB(Nyrs+ages-1);
NumericVector SB40_vecA(Nyrs);
NumericVector SB40_vecB(Nyrs);
NumericMatrix allelesA;
NumericMatrix allelesB;
NumericMatrix arichAmat;
NumericMatrix arichBmat;
NumericMatrix richnessA;
NumericMatrix richnessB;
NumericVector fst_vec(ngen+Nyrs);
NumericVector fst_vec2(ngen+Nyrs);
NumericVector Wrights_fst_vec(ngen+Nyrs);
NumericVector Wrights_fst_vec1(ngen+Nyrs);
NumericVector Wrights_fst_vec2(ngen+Nyrs);
NumericVector Wrights_simple(ngen+Nyrs);
NumericVector Wrights_simple1(ngen+Nyrs);
NumericVector Wrights_fst_vec_sub(ngen+Nyrs);
NumericVector Wrights_fst_vec_sub1(ngen+Nyrs);
NumericVector Wrights_simple_sub(ngen+Nyrs);

NumericVector arichsampA_vec(ngen+Nyrs);
NumericVector arichsampB_vec(ngen+Nyrs);
NumericVector arichsampA_ten(nSNP);
NumericVector arichsampB_ten(nSNP);
NumericVector ninecomp(9);
NumericVector diffsamp(10);
NumericVector nmig_vec(Nyrs);
NumericVector ps;
NumericVector s2;
NumericMatrix arichAmat_sub620;
NumericMatrix arichBmat_sub620;
NumericMatrix arichAmat_sub640;
NumericMatrix arichBmat_sub640;
NumericMatrix arichAmat_sub680;
NumericMatrix arichBmat_sub680;
//NumericVector arichAmat_sub;
//NumericVector arichBmat_sub;
NumericMatrix talmat2;
NumericMatrix talmat_sub2;
NumericVector NtotalA2;
NumericVector count_sub(nSNP);//list of number of unique alleles at each locus


if (pop1_prop==1){ninecomp=NumericVector::create(1,1,1,1,1,1,1,1,1);
if(samptype==1){diffsamp=rep(nsamp,10);}else{
diffsamp=NumericVector::create(200,180,155,130,110,90,65,40,20,10);}}
if (pop1_prop==0.2){ninecomp=NumericVector::create(1,2,2,2,2,2,2,2,2);
if(samptype==1){diffsamp=rep(nsamp,10);}else{
diffsamp=NumericVector::create(324,177,118,88,71,59,50,44,39,29);}}
if (pop1_prop==0.5){ninecomp=NumericVector::create(1,2,1,1,2,2,2,2,2);
if(samptype==1){diffsamp=rep(nsamp,10);}else{
diffsamp=NumericVector::create(170,130,100,65,33,170,130,100,65,33);}}

if (pop1_prop<=0.1){ninecomp=NumericVector::create(2,2,2,2,2,2,2,2,2);
if(samptype==1){diffsamp=rep(nsamp,10);}else{
diffsamp=NumericVector::create(324,177,118,88,71,59,50,44,39,29);}}

if (pop1_prop==0.01){ninecomp=NumericVector::create(2,2,2,2,2,2,2,2,2);
if(samptype==1){diffsamp=rep(nsamp,10);}else{
diffsamp=NumericVector::create(324,177,118,88,71,59,50,44,39,29);}}

if (pop1_prop==0.05){ninecomp=NumericVector::create(2,2,2,2,2,2,2,2,2);
if(samptype==1){diffsamp=rep(nsamp,10);}else{
diffsamp=NumericVector::create(324,177,118,88,71,59,50,44,39,29);}}


for (int zcounter=0;zcounter<nsims;zcounter ++){  //change back to 10
Rcout<<"zcounter"<<std::endl;
Rf_PrintValue(wrap(zcounter));
List data(nSNP+1);
NumericVector MUS(nSNP);
NumericVector SSBA_vec(Nyrs);
NumericVector SSBB_vec(Nyrs);

int num_wrong_pops;
int num_wrong_pops60;

NumericVector which_wrong_pops;  //vector of the split locations
NumericMatrix R_l0_mat;//this stores R_l0 estimates for each wrong pop for 100 yers
NumericMatrix epsilon_vecmat;
NumericMatrix TAC_vec;
NumericVector which_wrong_pops60;  //vector of the split locations
NumericMatrix R_l0_mat60;//this stores R_l0 estimates for each wrong pop for 100 yers
NumericMatrix epsilon_vecmat60;
NumericMatrix TAC_vec60;

if (popdyn==1){
// sim script - this is the part that sets up alleles at each locus
for (int m=0; m < nSNP; m++){  // this is for the 10 loci we simulate	

NumericVector alleles;
double prob;
double dir;
double mu;
int pick;
NumericVector alleles1(2*N);
NumericVector PICKS(2*N);
mu=0.1*as<double>(rbeta(1,0.6,1));//here scale down by 1/100 so max is 0.01 OR not
MUS[m]=mu;  //this takes mu back to 0.01* the rbeta  
alleles=rep(200,2*N);  //cool! this works in C++
for (int k=0; k < Yrs; k++){  //repeat mutation and random mating each year

//here is the mutate part

for (int i=0; i < 2*N; i++){  //repeat mutation for each allele
prob=as<double>(runif(1,0,1));
dir=as<double>(runif(1,0,1));		
if(prob<mu&dir<0.5){	
alleles[i]=alleles[i]+1;
}else if (prob<mu&dir>0.5){
alleles[i]=alleles[i]-1;
}			
}

//now the random mating part - first choose a random allele
for (int j=0; j < 2*N; j++){  //repeat random mating 2N times
prob=as<double>(runif(1,0,(2*N)-1));
pick=int(prob);
if((prob+0.5)>=(int(prob)+1)){	
pick=int(prob)+1;
}
//take that random allele and save it
PICKS[j]=pick;
alleles1[j]=alleles[pick];
}
alleles=alleles1;
}//end of Yrs
data[m]=alleles1;
}//end of m loop
data[nSNP]=MUS;

// end of simRcpp script
}//end popdyn==1
// initial numbers section (if there was one recruit)

N_til_init_a[0]=1;
for (int i=1; i < (ages-1); i++){  //fill in N_til_init;
N_til_init_a[i]=N_til_init_a[i-1]*exp(-M2[i-1]);}
N_til_init_a[ages-1]=N_til_init_a[ages-2]*exp(-M2[ages-2])/(1-exp(-M2[ages-1]));

//establish the SBPR with 1 recruit and F40%
N_til_init_a40[0]=1;
for (int i=1; i < (ages-1); i++){  //fill in N_til_init;
N_til_init_a40[i]=N_til_init_a40[i-1]*exp(-(M2[i-1]+fishsel[i-1]*fishmort));}
N_til_init_a40[ages-1]=N_til_init_a40[ages-2]*exp(-(M2[ages-2]+fishsel[ages-2]*fishmort))/(1-exp(-(M2[ages-1]+fishsel[ages-1]*fishmort)));

//for (int i=1; i < (ages-1); i++){  //fill in N_til_init;
//N_til_init_a40[i]=N_til_init_a40[i-1]*Surv[i-1];}
//N_til_init_a40[ages-1]=N_til_init_a40[ages-2]*Surv[ages-2]+N_til_init_a40[ages-1]*Surv[ages-1];


// establish length, weight of males and females
for (int i=0; i < ages; i++){ 
Length[i]=Linf*(1-exp(-K*(age[i])));
Weight_b[i]=psi*pow((Linf*(1-exp(-K*(age[i])))),theta);
}  //0.001 is for kg to mt conversion

Weight_a[0]=0;Weight_a[1]=0.5508;Weight_a[2]=1.1016;Weight_a[3]=1.6524;
Weight_a[4]=2.2032;Weight_a[5]=2.7540;Weight_a[6]=3.3048;

// determine initial number of recruits and the SSB in each population (pops 1 and 2)
S_l0[0]=SSBinitA;  //S_l0 initial spawning stock biomass for pop1
S_l0[1]=SSBinitB;    //S_l0 initial spawning stock biomass pop2     
NumericVector tempvec1a = Weight_a*N_til_init_a*Q_a;
NumericVector tempvec1b = Weight_b*N_til_init_a*Q_b;

R_l0A[0]=2.0*S_l0[0]/std::accumulate(tempvec1a.begin(), tempvec1a.end(), 0.0); //initial number recruits pop1  NOTE: if only one pop, all are in popA
R_l0B[0]=2.0*S_l0[1]/std::accumulate(tempvec1b.begin(), tempvec1b.end(), 0.0); //init recruits pop2
double R_l0A_init=R_l0A[0]; //record the first recruits for beginning of mgmt_counter loop
double R_l0B_init=R_l0B[0];

//n_init_a1 is the number of fish in each age class
//this loop will also round to an even number in each age class so you can split into males and females
for (int i=0; i<ages; i++){
N_init_a1[i]=N_til_init_a[i]*R_l0A[0];
temp=int(N_init_a1[i]);
if((N_init_a1[i]+0.5)>=(int(N_init_a1[i])+1)){	
temp=int(N_init_a1[i])+1;}
N_init_a1[i]=temp;
if (temp%2!=0){			
temp2=as<double>(runif(1,-1,1));
if(temp2>=0){N_init_a1[i]=N_init_a1[i]+1;}
else{N_init_a1[i]=N_init_a1[i]-1;}
}				
N_init_a2[i]=N_til_init_a[i]*R_l0B[0];
temp=int(N_init_a2[i]);
if((N_init_a2[i]+0.5)>=(int(N_init_a2[i])+1)){	
temp=int(N_init_a2[i])+1;}
N_init_a2[i]=temp;
if (temp%2!=0){	
temp2=as<double>(runif(1,-1,1));
if(temp2>=0){N_init_a2[i]=N_init_a2[i]+1;}
else{N_init_a2[i]=N_init_a2[i]-1;}
}

}//end this for loop

N_init_a1_FIRST = clone(N_init_a1);
N_init_a2_FIRST = clone(N_init_a2);
N_init_a = N_init_a1+N_init_a2;

NumericVector tempvec_a = Weight_a*N_init_a*Q_a;
NumericVector tempvec_b = Weight_b*N_init_a*Q_b;
S_lyINIT_a=0.5*std::accumulate(tempvec_a.begin(), tempvec_a.end(), 0.0);//Initial Spawning bimass
S_lyINIT_b=0.5*std::accumulate(tempvec_b.begin(), tempvec_b.end(), 0.0);//Initial Spawning bimass
SPR_init_a=S_lyINIT_a/N_init_a[0]; //initial SPBR
SPR_init_b=S_lyINIT_b/N_init_a[0]; //initial SPBR
tempvec_a = Weight_a*N_til_init_a40*Q_a;
tempvec_b = Weight_b*N_til_init_a40*Q_b;
SPR_init_F40_a=0.5*std::accumulate(tempvec_a.begin(), tempvec_a.end(), 0.0);//Initial Spawning bimass
SPR_init_F40_b=0.5*std::accumulate(tempvec_b.begin(), tempvec_b.end(), 0.0);//Initial Spawning bimass

NumericVector temp6=(N_init_a1*Weight_a*Q_a);
double S_lyINITA=0.5*std::accumulate(temp6.begin(),temp6.end(),0.0);
NumericVector temp7=(N_init_a2*Weight_b*Q_b);
double S_lyINITB=0.5*std::accumulate(temp7.begin(),temp7.end(),0.0);
double SPR_initA;  //initial SBPR popA
double SPR_initB;  //initial SBPR popB

pop1_size = std::accumulate(N_init_a1.begin(),N_init_a1.end(),0.0);
pop2_size = std::accumulate(N_init_a2.begin(),N_init_a2.end(),0.0);

int N_init=pop1_size+pop2_size;

if (popdyn==1){
// evolve script - takes alleles for all individuals, splits it (or not in the case of one population) and then evolves population(s) for ngen generations
NumericVector transferdata(N_init);

//size of pop1 and pop2 (out of the original size)
NumericVector alleles_pop1_before(2*pop1_size);
NumericVector alleles_pop2_before(2*pop2_size);
NumericVector alleles_pop1_after(2*pop1_size); //after random mating
NumericVector alleles_pop2_after(2*pop2_size); //after random mating
NumericVector PICKS_pop1(2*pop1_size);
NumericVector PICKS_pop2(2*pop2_size);
//here we extract elements of the data list (with 10 microsatellites)
NumericVector selectA(2*pop1_size);
NumericVector selectB(2*pop2_size);
int pick3;
for (int i=0; i < 2*pop1_size; i++){  //create a list of random numbers 0-1999 twice as big as pop1_size to randomly select alleles from sim
prob=as<double>(runif(1,0,(2*N)-1)); //works, I checked it
pick3=int(prob);
if((prob+0.5)>=(int(prob)+1)){	
pick3=int(prob)+1;
}
selectA[i]=pick3;}

for (int i=0; i < 2*pop2_size; i++){  //create a list of random numbers 0-1999 twice as big as pop2_size to randomly select alleles from sim
prob=as<double>(runif(1,0,(2*N)-1));
pick3=int(prob);
if((prob+0.5)>=(int(prob)+1)){	
pick3=int(prob)+1;
}
selectB[i]=pick3;}

for (int i=0; i < nSNP; i++){  
transferdata = data[i];
NumericVector tx1(2*pop1_size);
NumericVector tx2(2*pop2_size);
for (int j=0; j < 2*pop1_size; j++){
tx1[j]=transferdata[selectA[j]];
}
for (int j=0; j < 2*pop2_size; j++){
tx2[j]=transferdata[selectB[j]];
}
pop1[i]=tx1;
pop2[i]=tx2;
}

for (int k=0; k < ngen; k++){ //for ngen generations
for (int j=0; j<nSNP; j++){ //for all 10 microsatellites
//mutate alleles
alleles_pop1_before = pop1[j]; 
for (int i=0; i < 2*pop1_size; i++){  //pop1 mutation for each allele
prob=as<double>(runif(1,0,1));  
dir=as<double>(runif(1,0,1));		
if(prob<(MUS[j])&dir<0.5){	
alleles_pop1_before[i]=alleles_pop1_before[i]+1;
}else if (prob<(MUS[j])&dir>0.5){
alleles_pop1_before[i]=alleles_pop1_before[i]-1;
}			
}	

alleles_pop2_before = pop2[j];		 
for (int i=0; i < 2*pop2_size; i++){  //pop2 mutation for each allele
prob=as<double>(runif(1,0,1));  
dir=as<double>(runif(1,0,1));		
if(prob<(MUS[j])&dir<0.5){	
alleles_pop2_before[i]=alleles_pop2_before[i]+1;
}else if (prob<(MUS[j])&dir>0.5){
alleles_pop2_before[i]=alleles_pop2_before[i]-1;
}			
}

//random mating 
for (int i=0; i < 2*pop1_size; i++){  //random mating for each allele
prob=as<double>(runif(1,-0.5,(2*(rmate*pop1_size)-0.5))); 
pick1=int(prob);
if((prob+0.5)>=(int(prob)+1)){	
pick1=int(prob)+1;
}
PICKS_pop1[i]=pick1;  //also seems to work, values from zero to (2*pop1_size)-1
alleles_pop1_after[i]=alleles_pop1_before[pick1];
}
//PICKS_pop1 is a set of random alleles from 0 to (2*pop1_size)-1 who will reproduce
for (int i=0; i < 2*pop2_size; i++){  //random mating for each allele
prob=as<double>(runif(1,-0.5,(2*(rmate*pop2_size)-0.5)));  //was 2*pop2_size-0.5
pick2=int(prob);
if((prob+0.5)>=(int(prob)+1)){	
pick2=int(prob)+1;
}
PICKS_pop2[i]=pick2;
alleles_pop2_after[i]=alleles_pop2_before[pick2];
}

//for each of j microsatellites, I have mutated them, and then selected some
//to reproduce next year. They are now stored in alleles_pop1_after and alleles_pop2_after
//migration  removed for now in order to achieve high enough Fst values between populations
//int nmigevol;

//nmigevol=nmig;
//if(nmigevol>0){
//for (int m=0; m<nmigevol; m++){// for each migrant
//mover1=alleles_pop1_after(PICKS_pop1[m]);
//alleles_pop1_after.erase(PICKS_pop1[m]);
//mover2=alleles_pop2_after(PICKS_pop2[m]);
//alleles_pop2_after.erase(PICKS_pop2[m]);
//alleles_pop1_after.insert(PICKS_pop1[m],mover2);
//alleles_pop2_after.insert(PICKS_pop2[m],mover1);
//} // for each migrant
//}  //if nmig>0

NumericVector clone_alleles_pop1_after=clone(alleles_pop1_after);
NumericVector clone_alleles_pop2_after=clone(alleles_pop2_after);
pop1[j]=clone_alleles_pop1_after;  //update the list of alleles
pop2[j]=clone_alleles_pop2_after;
}  //for each of 10 microsatellites

Rcout<<"k in evolve"<<std::endl;
Rf_PrintValue(wrap(k));

genotypes[0]=pop1;
genotypes[1]=pop2;
//end of evolve


//**********************************************************************************************
//TRY FST CALCS HERE

if (k%200==0){
allelesA = NumericMatrix((2*pop1_size),nSNP); //is it really 2*? yes
allelesB = NumericMatrix((2*pop2_size),nSNP); //is it really 2*? yes
arichAmat=NumericMatrix(2*pop1_size,nSNP);
arichBmat=NumericMatrix(2*pop2_size,nSNP);
richnessA=NumericMatrix(Nyrs,nSNP); //holds each year allelic richness for each microsat
richnessB=NumericMatrix(Nyrs,nSNP);
// now fill up allelesA with genotypes (one in each column)

List temp_genA;
List temp_genB;
NumericVector temp3;

//CHECK WHETERH THIS WORKS FOR 1 population
temp_genA=genotypes[0];
for (int i=0; i<nSNP; i++){
temp3=temp_genA[i];
allelesA(_,i)=temp3;  //each column of allelesA has all alleles for that usat
}
temp_genB=genotypes[1];
for (int i=0; i<nSNP; i++){
temp3=temp_genB[i];
allelesB(_,i)=temp3;
}
//make a matrix to hold genotypes for keeping track of allelic richness
//first year popA and popB

arichAmat=allelesA;  //arichAmat is all alleles males and females from population first year
arichBmat=allelesB;  //use allelesA same thing I think

//***********************************************************************************************
//Fst (for the entire population)  if nsamp is bigger than pop size, will give warning
//create a matrix with unique alleles for each microsatellite in each column

NumericMatrix tallyBOTHmat(2*pop2_size+2*pop1_size,nSNP);
int some_integerBOTH = 0;
NumericMatrix arichBOTHmat(2*pop2_size+2*pop1_size,nSNP);//all alleles both pops year 1

for (int i=0;i<2*pop1_size;i++){
arichBOTHmat(i,_)=arichAmat(i,_);
}
for (int i=0;i<2*pop2_size;i++){
arichBOTHmat(i+2*pop1_size,_)=arichBmat(i,_);
}

//get unique alleles for both populations
for (int rich = 0; rich<nSNP;rich++){
NumericMatrix::Column inputBOTH = arichBOTHmat(_,rich); 
NumericVector xBOTH = clone<NumericVector>(inputBOTH);
NumericVector tally_allelesBOTH(2*pop2_size+2*pop1_size); //this is a vector that holds the unique alleles
int tallycounter=1;
//*******
int nallsBOTH=0;
typedef std::map<double,int> imap ;
typedef imap::value_type pair ;
imap index ;
int n = xBOTH.size() ;
double current, previous = xBOTH[0] ;
index.insert( pair( previous, 0 ) );

imap::iterator it = index.begin() ;
for( int i=1; i<n; i++){
current = xBOTH[i] ;
if( current == previous ){
xBOTH[i] = current + ( ++(it->second) / 100.0 ) ;
} else {
it = index.find(current) ;
if( it == index.end() ){
it = index.insert(
current > previous ? it : index.begin(),
pair( current, 0 )
) ;
} else {
xBOTH[i] = current + ( ++(it->second) / 100.0 ) ;
}
previous = current ;
}

if (xBOTH[i]-inputBOTH[i]==0)
{nallsBOTH=nallsBOTH+1;
tally_allelesBOTH[tallycounter]=inputBOTH[i];
tallycounter=tallycounter+1;
}
}
tally_allelesBOTH[0]=inputBOTH[0];

//*****************************
tallyBOTHmat(_,rich)=tally_allelesBOTH;
} //end of rich BOTH loop (over 10 microsats)

List tallist(nSNP); //list of 10 that will hold allele frequencies for each of 10 microsatellites.
//First column is the actual allele size, second is freq from popA, third is freq popB
NumericVector count(nSNP);//list of number of unique alleles at each locus
//now set up a table to get allele frequencies for popA and popB
for (int rich=0;rich<nSNP;rich++){
NumericMatrix::Column inputA=arichAmat(_,rich);
NumericMatrix::Column inputB=arichBmat(_,rich);
NumericMatrix::Column tal=tallyBOTHmat(_,rich); 

for(int i=0;i<tal.size();i++){  //figure out how many unique alleles are at this microsat (from tal)
if(tal[i]!=0){count[rich]=count[rich]+1;}
}

NumericVector Counter=clone(count);
int counter=Counter[rich];
NumericVector taltrunc(counter);  //will hold all unique alleles from this microsat (both pops)
NumericMatrix tallyAB(counter,3);  //matrix that has alleles, freq at popA, freq at popB
NumericVector howmanyA(counter);  //number of alleles for this population at each allele
NumericVector howmanyB(counter);

for(int i=0;i<counter;i++){  //counter is the number of unique alleles at a locus
taltrunc[i]=tal[i];

int counterA=0;  //a counter for number of unique alleles at each locus
int counterB=0;

for (int j=0;j<2*pop1_size;j++){
if (inputA[j]==taltrunc[i])//go through all alleles to see how many match this unique one
{counterA=counterA+1;}
}
howmanyA[i]=counterA;

for (int j=0;j<2*pop2_size;j++){
if (inputB[j]==taltrunc[i])
{counterB=counterB+1;}
}

howmanyB[i]=counterB;

} //end of counter

tallyAB(_,0)=taltrunc;
tallyAB(_,1)=howmanyA/(2*pop1_size);
tallyAB(_,2)=howmanyB/(2*pop2_size);

tallist[rich]=tallyAB;

}//end of rich

//create talmat, which has unique alleles first column then freqs for popA and B
NumericMatrix talmat(std::accumulate(count.begin(),count.end(),0.0),5);

for (int i=0;i<nSNP;i++){
int talcount=0;
NumericMatrix taltmp(count[i],3);
taltmp=as<SEXP>(tallist[i]);
for (int j=std::accumulate(count.begin(),count.begin()+i,0.0);j<std::accumulate(count.begin(),count.begin()+i+1,0.0);j++)
{
talmat(j,_)=taltmp(talcount,_);
talcount=talcount+1;
}
}
//aha! talmat is a matrix with all loci. Genious!


//talmat is the raw material to find Fst!!!
//return(wrap(talmat));
//GET FST this function finds Fst between 2 populations, multiple loci and multiple alleles.
double n_sampA=pop1_size;
double n_sampB=pop2_size;
double n_bar=0.5*(n_sampA+n_sampB); //sample size - can change this.
double r=2;
double C=0;
NumericVector p_bar(talmat.nrow());
NumericMatrix::Column sampmatA =talmat(_,1);
NumericMatrix::Column sampmatB =talmat(_,2); 
s2=NumericVector(talmat.nrow());
NumericVector h_bar(talmat.nrow());
NumericVector ones(talmat.nrow());
ones=rep(1,talmat.nrow());
p_bar=(n_sampA*sampmatA+n_sampB*sampmatB)/(2*n_bar);  //each entry is average sample frequency of an allele
s2=pow((sampmatB-p_bar),2)+pow((sampmatA-p_bar),2);//sample variance of allele freqs over pops (for each allele)

//calculate mean heterozygosity (He=1-sum(all allele frequencies squared)), aka remove the homozygotes
h_bar=((n_sampA*2*sampmatA*(ones-sampmatA))+(n_sampB*2*sampmatB*(ones-sampmatB)))/(2*n_bar);  //this takes the mean of 2 populations

//here heterozygosity assumes HWE.

double nc=((r*n_bar)-((pow(n_sampA,2)/(r*n_bar))+(pow(n_sampB,2)/(r*n_bar))))/(r-1);  //same as n_bar

NumericVector a=(n_bar/nc)*(s2-(1/(n_bar-1))*((p_bar*(1-p_bar))-((r-1)/r)*s2-(0.25*h_bar)));
NumericVector dvec=((2*n_bar)/((2*n_bar)-1))*((p_bar*(1-p_bar))-((r-1)/r)*s2);
NumericVector b=(n_bar/(n_bar-1))*((p_bar*(1-p_bar))-((r-1)/r)*s2-(h_bar*(2*n_bar-1)/(4*n_bar)));
NumericVector c=h_bar/2;

NumericVector aplusdvec=a+b+c;
double fst=std::accumulate(a.begin(),a.end(),0.0)/(std::accumulate(aplusdvec.begin(),aplusdvec.end(),0.0)); //the better one
fst_vec[k]=fst;
}//second k loop
}//if k%100
}//end popdyn==1

//END TRY FST CALCS
//**********************************************************************************************
//mgmt_counter 
//0=combined
//1=separated
//2=separated genetics test
//3=combined then separated
//4=tier 5
//5=combined then no fishing
//6=combined then tier5

//now take note of whether there are any false boundaries when the true boundary is not found
//here you are looking at the significant test results (mistakes)
//I am only interested in these if you do not correctly identify true split
//and only if there really are 2 populations
//but will calculate them here for all cases and then filter later

for (int mgmt_counter=0;mgmt_counter<mmt;mgmt_counter ++){
N_init_a1=N_init_a1_FIRST;
N_init_a2=N_init_a2_FIRST;

//genrule is the integer that tells you whether genetics test found the correct split site.

for (int k=0; k < Nyrs; k++){  //make k to be 110		******should be Nyrs but 1 for now
//calculate age structure next year

temp8=(N_init_a1*Weight_a*Q_a);
S_lyA=0.5*std::accumulate(temp8.begin(),temp8.end(),0.0);  //spawning biomass of popA

temp9=(N_init_a2*Weight_b*Q_b);
S_lyB=0.5*std::accumulate(temp9.begin(),temp9.end(),0.0);  //spawning biomass of popB

SSBA_vec[k]=S_lyA;
SSBB_vec[k]=S_lyB;

if(N_init_a1[0]==0){SPRAtrue[k]=0;}else{SPRAtrue[k]=S_lyA/N_init_a1[0];}//true spawning biomass per recruit in each year
if(N_init_a2[0]==0){SPRBtrue[k]=0;}else{SPRBtrue[k]=S_lyB/N_init_a2[0];}

if(k>=1){
eta_jyA=as<double>(rnorm(1,0,sqrt(sigma_ja2A)));
eta_jyB=as<double>(rnorm(1,0,sqrt(sigma_ja2B)));
eta_jy=as<double>(rnorm(1,0,sqrt(sigma_ja2)));
epsilonA[k]=rho*epsilonA[k-1]+pow((1-pow(rho,2)),0.5)*eta_jyA;
epsilonB[k]=rho*epsilonB[k-1]+pow((1-pow(rho,2)),0.5)*eta_jyB;
epsilon[k]=rho*epsilon[k-1]+pow((1-pow(rho,2)),0.5)*eta_jy;		
} //end of k>=1 loop

N_hat_aA = N_init_a1*exp(epsilonA[k]-sigma_ja2A/2);
N_hat_aB = N_init_a2*exp(epsilonB[k]-sigma_ja2B/2);
N_hat_a=(N_init_a1+N_init_a2)*exp(epsilon[k]-sigma_ja2/2);
psiA=as<double>(rnorm(1,0,sqrt(sigma_js2A)));
psiB=as<double>(rnorm(1,0,sqrt(sigma_js2B)));
psi_star=as<double>(rnorm(1,0,sqrt(sigma_js2)));

R_l0A_hat[0]=N_hat_aA[0];//same correction applied here as to total numbers
R_l0B_hat[0]=N_hat_aB[0];

if(N_hat_aA[0]==0){S_hat_syjA[k]=0;}else{
temp10= Weight_a*N_hat_aA*Q_a;
S_hat_syjA[k]=0.5*std::accumulate(temp10.begin(),temp10.end(),0.0)*exp(psiA-(sigma_js2A/2))/N_hat_aA[0];}  //simulated estimate of spawning biomass per recruit
SSB_hatA[k]= 0.5*std::accumulate(temp10.begin(),temp10.end(),0.0)*exp(psiA-(sigma_js2A/2));		    
if(N_hat_aB[0]==0){S_hat_syjB[k]=0;}else{ 
temp11 = Weight_b*N_hat_aB*Q_b;  
S_hat_syjB[k]=0.5*std::accumulate(temp11.begin(),temp11.end(),0.0)*exp(psiB-(sigma_js2B/2))/N_hat_aB[0];}
SSB_hatB[k]= 0.5*std::accumulate(temp11.begin(),temp11.end(),0.0)*exp(psiB-(sigma_js2B/2));        
if(N_hat_a[0]==0){S_hat_syj[k]=0;}else{
temp12 = Weight_a*N_hat_a*Q_a;
S_hat_syj[k]=0.5*std::accumulate(temp12.begin(),temp12.end(),0.0)*exp(psi-(sigma_js2/2))/N_hat_a[0];}
SSB_hat[k]= 0.5*std::accumulate(temp12.begin(),temp12.end(),0.0)*exp(psi-(sigma_js2/2));        

SPR_initA=SPR_init_F40_a*(std::accumulate(R_l0A_hat.begin(),R_l0A_hat.begin()+k,0.0)/k);//maybe get rid of this and other bad measures using it
SPR_initB=SPR_init_F40_b*(std::accumulate(R_l0B_hat.begin(),R_l0B_hat.begin()+k,0.0)/k);
SB40_vecA[k]=SPR_init_F40_a*(std::accumulate(R_l0A_hat.begin(),R_l0A_hat.begin()+k,0.0)/k);
SB40_vecB[k]=SPR_init_F40_b*(std::accumulate(R_l0B_hat.begin(),R_l0B_hat.begin()+k,0.0)/k);


FishA=as<double>(rlnorm(1,-0.39,0.46));//this if FullF for B20%  
FishB=as<double>(rlnorm(1,-0.39,0.46));//this is FullF for B40%


if (k<500){//this is just mean fishsel
fishselB[0]=0;fishselB[1]=0.818;fishselB[2]=0.974;fishselB[3]=0.895;fishselB[4]=0.847;fishselB[5]=0.848;fishselB[6]=0.848;
}

if (k>=500){//this is low fishsel
fishselB[0]=0;fishselB[1]=0.4624;fishselB[2]=0.6072;fishselB[3]=0.5354;fishselB[4]=0.5110;fishselB[5]=0.5282;fishselB[6]=0.5282;
}


FishAvec[k]=FishA;
FishBvec[k]=FishB;

temp15=Weight_a*N_init_a1*((fishsel*FishA)/(fishsel*FishA+M2))*(1-exp(-(fishsel*FishA+M2)));
TACA[k]=std::accumulate(temp15.begin(),temp15.end(),0.0);
temp16=Weight_b*N_init_a2*((fishselB*FishB)/(fishselB*FishB+M2))*(1-exp(-(fishselB*FishB+M2)));
TACB[k]=std::accumulate(temp16.begin(),temp16.end(),0.0);

//get spawning fish numbers
spawnersA = Q_a*N_init_a1*0.5;  
spawnersB = Q_b*N_init_a2*0.5;

for (int i=0; i<ages; i++){
temp=int(spawnersA[i]);
if((spawnersA[i]+0.5)>=(int(spawnersA[i])+1)){	
temp=int(spawnersA[i])+1;
}
spawnersA_rounded[i]=temp;
}
for (int i=0; i<ages; i++){
temp=int(spawnersB[i]);
if((spawnersB[i]+0.5)>=(int(spawnersB[i])+1)){	
temp=int(spawnersB[i])+1;
}
spawnersB_rounded[i]=temp;
}

//make sure there are not more spawners than actual fish.
for(int i=0;i<ages;i++){
if((N_init_a1[i]/2)<spawnersA_rounded[i]){
spawnersA_rounded[i]=N_init_a1[i]/2;
}
if((N_init_a2[i]/2)<spawnersB_rounded[i]){
spawnersB_rounded[i]=N_init_a2[i]/2;
}
}

NspawnersA[k]=std::accumulate(spawnersA_rounded.begin(),spawnersA_rounded.end(),0.0); //males and females same number of spawners as in here
NspawnersB[k]=std::accumulate(spawnersB_rounded.begin(),spawnersB_rounded.end(),0.0);

//generate delta which is the variance for recruitment
double sigma_jr2_tot=log(pow(sigma_jr2,2)+1);
double sigma_rec=as<double>(rnorm(1,0,sqrt(sigma_jr2_tot)));

//N_aA and N_aB are numbers next year*****could double check this part  
N_aA[0]=((4*h*R_l0A[0]*S_lyA)/(S_l0[0]*(1-h)+S_lyA*(5*h-1)))*exp(sigma_rec-(sigma_jr2/2));  //S_l0[0] are spawners in popA //BH spawner recruit relationship
N_aB[0]=((4*h*R_l0B[0]*S_lyB)/(S_l0[1]*(1-h)+S_lyB*(5*h-1)))*exp(sigma_rec-(sigma_jr2/2));  //S_l0[1] are spawners in popB

Rcout<<"first migration"<<std::endl;
//HERE IS MIGRATION IN POP DYN MODEL
//if (k%4==0){
//for (int i=0;i<nmig;i++){
//N_aA[0]=N_aA[0]+2;
//}
//}//end of if k%=4==0

if(NspawnersA[k]==0){N_aA[0]=0;}
if(NspawnersB[k]==0){N_aB[0]=0;}

for(int i=1; i<(ages-1); i++){
N_aA[i]=N_init_a1[i-1]*exp(-(fishsel[i-1]*FishA+M2[i-1]));  //correct because you want the selectivity of the next youngest age
N_aB[i]=N_init_a2[i-1]*exp(-(fishselB[i-1]*FishB+M2[i-1]));  //and fishsel now goes from 21 to 1
//N_aA[i]=N_init_a1[i-1]*Surv[i-1];  //correct because you want the selectivity of the next youngest age
//N_aB[i]=N_init_a2[i-1]*Surv[i-1];  //and fishsel now goes from 21 to 1
}
N_aA[ages-1]=N_init_a1[ages-2]*exp(-(fishsel[ages-2]*FishA+M2[ages-2]))+N_init_a1[ages-1]*exp(-(fishsel[ages-1]*FishA+M2[ages-1]));   //dont forget about plus group
N_aB[ages-1]=N_init_a2[ages-2]*exp(-(fishselB[ages-2]*FishB+M2[ages-2]))+N_init_a2[ages-1]*exp(-(fishselB[ages-1]*FishB+M2[ages-1]));

//N_aA[ages-1]=N_init_a1[ages-2]*Surv[ages-2]+N_init_a1[ages-1]*Surv[ages-1];   //dont forget about plus group
//N_aB[ages-1]=N_init_a2[ages-2]*Surv[ages-2]+N_init_a2[ages-1]*Surv[ages-1];


//now round N_aA and N_aB to integers
for (int i=0; i<ages; i++){
temp18=int(N_aA[i]);
if((N_aA[i]+0.5)>=(int(N_aA[i])+1)){	
temp18=int(N_aA[i])+1;
}
N_aA[i]=temp18;
if (temp18%2!=0){		
temp2=as<double>(runif(1,-1,1));
if(temp2>=0){N_aA[i]=N_aA[i]+1;}
else{N_aA[i]=N_aA[i]-1;
}
}
}
for (int i=0; i<ages; i++){
temp19=int(N_aB[i]);
if((N_aB[i]+0.5)>=(int(N_aB[i])+1)){	
temp19=int(N_aB[i])+1;}
N_aB[i]=temp19;
if (temp19%2!=0){		
temp2=as<double>(runif(1,-1,1));
if(temp2>=0){N_aB[i]=N_aB[i]+1;}
else{N_aB[i]=N_aB[i]-1;
}
}
}


//how many individuals are there? 
R_l0A[k]=N_init_a1[0];
R_l0B[k]=N_init_a2[0];
R_l0A_hat[k]=N_hat_aA[0]; //incorporate error into recruitment observation
R_l0B_hat[k]=N_hat_aB[0];

NtotalA[k] = std::accumulate(N_init_a1.begin(),N_init_a1.end(),0.0);
NtotalB[k] = std::accumulate(N_init_a2.begin(),N_init_a2.end(),0.0);




allA.row(k)=rev(N_init_a1);  //each row is the number of each age during that year  - reversed so it is 21-1 for next section
allB.row(k)=rev(N_init_a2); 
allspawnA.row(k)=rev(spawnersA_rounded);  //each row is the number of each age during that year - reversed so it is 21-1 for next section
allspawnB.row(k)=rev(spawnersB_rounded);

//put this at the end 
N_init_a1=clone(N_aA);
N_init_a2=clone(N_aB);

//create separate matrices for each of the mgmt cases. Can return for separated case (left out here).

out3(zcounter,k)=R_l0A[k];//true recruitment
out4(zcounter,k)=R_l0B[k];
out21(zcounter,k)=0.4*SSBA_vec[k]/SB40_vecA[k];//true SSB div by our GThompson est of B40
out22(zcounter,k)=0.4*SSBB_vec[k]/SB40_vecB[k];
out23(zcounter,k)=SSBA_vec[k];//just ssb
out24(zcounter,k)=SSBB_vec[k];


}  //end of k<Nyrs loop

NtotalA2=clone(NtotalA);

//try calculating Ne here work from 1 to <ages> here
NumericVector x_age(ages);
NumericVector S_xA(ages);
NumericVector S_xB(ages);
NumericVector bxA(ages);
NumericVector bxB(ages);
NumericVector lxA(ages);
NumericVector lxB(ages);
NumericVector b_primexA(ages);
NumericVector b_primexB(ages);
NumericVector bxlxA(ages);
NumericVector bxlxB(ages);
NumericVector BxA(ages);
NumericVector BxB(ages);
double L_popA; //generation length
double L_popB; //generation length
NumericVector L_popAvec(ages);
NumericVector L_popBvec(ages);
NumericVector L_popA_vec(Nyrs);
NumericVector L_popB_vec(Nyrs);
NumericVector VkA_vec(Nyrs);
NumericVector VkB_vec(Nyrs);
NumericVector k_barA(ages);
NumericVector k_barB(ages);
NumericVector VxA(ages);
NumericVector VxB(ages);
NumericVector DxA(ages);
NumericVector DxB(ages);
NumericVector k_barADxA(ages);
NumericVector k_barBDxB(ages);
NumericVector SSDIxA(ages);
NumericVector SSDIxB(ages);
NumericVector delxA(ages);
NumericVector delxB(ages);
NumericVector SSDGxA(ages);
NumericVector SSDGxB(ages);
NumericVector SSDxA(ages);
NumericVector SSDxB(ages);
NumericVector newnumsA(ages);//calculate new numbers based on initial nums and lx.
NumericVector newnumsB(ages);
double k_bardotA;
double k_bardotB;
double VkA;
double VkB;

for (int k=0;k<Nyrs;k++){  //put in Nyrs
x_age=age+1;
N_init_a1_star = rev(allA.row(k));  //0 yr olds to <ages> yr olds here
N_init_a2_star = rev(allB.row(k));
spawnersA_rounded_star = rev(allspawnA.row(k));
spawnersB_rounded_star = rev(allspawnB.row(k));
FishA = FishAvec[k];
FishB = FishBvec[k];
//S_xA = 1-(fishsel*FishA+M);
//S_xB = 1-(fishsel*FishB+M);
S_xA=rep((1-meanFullF),ages);
S_xA[(ages-1)]=0;
S_xB=rep((1-meanFullF),ages);
S_xB[(ages-1)]=0;
for (int i=0;i<ages;i++){  //make sure they are positive (neg can happen when fishing gets too high)
if (S_xA[i]<0){S_xA[i]=0;}
if (S_xB[i]<0){S_xB[i]=0;}
}
lxA[0]=1;
lxB[0]=1;
for (int i=1;i<ages;i++){
lxA[i] = S_xA[i-1]*lxA[i-1];
lxB[i] = S_xB[i-1]*lxB[i-1];
}
newnumsA=0.5*N_init_a1_star[0]*lxA;
newnumsB=0.5*N_init_a2_star[0]*lxB;
for(int i=0;i<ages;i++){
if (NspawnersA[k]*newnumsA[i]>0)
{bxA[i] = spawnersA_rounded_star[i]*N_init_a1_star[0]/(2*NspawnersA[k]*newnumsA[i]);}//div by2 because want half the number of N1s here
else{bxA[i]=0;}
if (NspawnersB[k]*newnumsB[i]>0)
{bxB[i] = spawnersB_rounded_star[i]*N_init_a2_star[0]/(2*NspawnersB[k]*newnumsB[i]);}
else{bxB[i]=0;}
}
bxlxA=bxA*lxA;
bxlxB=bxB*lxB;
for(int i=0;i<ages;i++){
if(std::accumulate(bxlxA.begin(),bxlxA.end(),0.0)>0){
b_primexA[i] = 2*bxA[i]/std::accumulate(bxlxA.begin(),bxlxA.end(),0.0);}
else{b_primexA[i] =0;}
if(std::accumulate(bxlxB.begin(),bxlxB.end(),0.0)>0){
b_primexB[i] = 2*bxB[i]/std::accumulate(bxlxB.begin(),bxlxB.end(),0.0);}
else{b_primexB[i]=0;}
}
BxA=b_primexA*newnumsA;
BxB=b_primexB*newnumsB;
L_popAvec = x_age*BxA;
L_popBvec = x_age*BxB;
if (N_init_a1_star[0]>0){
L_popA = std::accumulate(L_popAvec.begin(),L_popAvec.end(),0.0)/N_init_a1_star[0];}
else {L_popA=0;}
if (N_init_a2_star[0]>0){
L_popB = std::accumulate(L_popBvec.begin(),L_popBvec.end(),0.0)/N_init_a2_star[0];}
else{L_popB=0;}

L_popA_vec[k]=L_popA;
L_popB_vec[k]=L_popB;

k_barA[0]=0;
k_barB[0]=0;
for (int i=1;i<ages;i++){
k_barA[i] = b_primexA[i]+k_barA[i-1];
k_barB[i] = b_primexB[i]+k_barB[i-1];
}

VxA=k_barA;//just for males or females (but they have the same value) assume poisson variance
VxB=k_barB;


for (int i=0;i<(ages-1);i++){

if (newnumsA[i]-newnumsA[i+1]>=0){DxA[i] =newnumsA[i]-newnumsA[i+1];}
else{DxA[i]=0;}
if (newnumsB[i]-newnumsB[i+1]>=0){DxB[i] =newnumsB[i]-newnumsB[i+1];}
else{DxB[i]=0;}
}
k_barADxA=k_barA*DxA;
k_barBDxB=k_barB*DxB;
if (N_init_a1_star[0]>0){
k_bardotA = std::accumulate(k_barADxA.begin(),k_barADxA.end(),0.0)/N_init_a1_star[0];}
else{k_bardotA=0;}//divide by number of newborns
if (N_init_a2_star[0]>0){
k_bardotB = std::accumulate(k_barBDxB.begin(),k_barBDxB.end(),0.0)/N_init_a2_star[0];}
else{k_bardotB=0;}
delxA=k_barA-k_bardotA;
delxB=k_barB-k_bardotB;

SSDIxA=DxA*VxA;
SSDIxB=DxB*VxB;

SSDGxA=DxA*delxA*delxA;
SSDGxB=DxB*delxB*delxB;

SSDxA=SSDIxA+SSDGxA;
SSDxB=SSDIxB+SSDGxB;
double SSDtA = std::accumulate(SSDxA.begin(),SSDxA.end(),0.0);
double SSDtB = std::accumulate(SSDxB.begin(),SSDxB.end(),0.0);
if (N_init_a1_star[0]>0){VkA = SSDtA/N_init_a1_star[0];} //just for females or males but the values are the same
else{VkA=0;}
if (N_init_a2_star[0]>0){VkB = SSDtB/N_init_a2_star[0];}
else{VkB=0;}
//VkA and VkB are total variance 
VkA_vec[k]=VkA;
VkB_vec[k]=VkB;


double NeA = (4*N_init_a1_star[0]*L_popA)/(VkA+2);
double NeB = (4*N_init_a2_star[0]*L_popB)/(VkB+2);

//if(NtotalA[k]>1){
//NeAvec[k]=NeA/NtotalA[k];}
//else{NeAvec[k]=0;}
//if(NtotalB[k]>1){
//NeBvec[k]=NeB/NtotalB[k];}
//else{NeBvec[k]=0;}

//calculate vkdot from males and females in pop and also teh VkA an dVkB.
//then find Ne.
NeAvec[k]=NeA;
NeBvec[k]=NeB;

}//end of k loop for finding Ne.


out7(zcounter,_) = TACA;
out8(zcounter,_)= TACB;
out9(zcounter,_)=FishAvec;
out10(zcounter,_)=FishBvec;

out14(zcounter,_)=TAC;
out15(zcounter,_)=NtotalA;
out16(zcounter,_)=NtotalB;

if (popdyn==1){
//Set up lists which will hold genotypes for males, females each pop
//List popA(Nyrs+ages);
//List popB(Nyrs+ages);

//popA and B lists of 131 matrices. First 21 will be length of N_init_a1_FIRST (age 21 to 1), the rest will be number of offspring for that year
NumericMatrix matsliceA;  //popA males
NumericMatrix matsliceB;  //popA females
//popA set up alleles for each individual from allelesA
for(int i=0;i<ages;i++){
matsliceA=NumericMatrix(N_init_a1_FIRST[i],nSNP);
matsliceB=NumericMatrix(N_init_a1_FIRST[i],nSNP);
for (int j=std::accumulate(N_init_a1_FIRST.begin(),N_init_a1_FIRST.begin()+i,0.0);
j<std::accumulate(N_init_a1_FIRST.begin(),N_init_a1_FIRST.begin()+i+1,0.0);j++){
matsliceA(j-std::accumulate(N_init_a1_FIRST.begin(),N_init_a1_FIRST.begin()+i,0.0),_)=allelesA(j,_);}
for (int j=std::accumulate(N_init_a1_FIRST.begin(),N_init_a1_FIRST.begin()+i,0.0)+pop1_size;
j<std::accumulate(N_init_a1_FIRST.begin(),N_init_a1_FIRST.begin()+i+1,0.0)+pop1_size;j++){
matsliceB(j-std::accumulate(N_init_a1_FIRST.begin(),N_init_a1_FIRST.begin()+i,0.0)-pop1_size,_)=allelesA(j,_);}

List temp5(2); temp5[0]=matsliceA;temp5[1]=matsliceB;
popA[ages-1-i]=temp5;
}

//pop2 set up alleles for each individual from allelesB 
NumericMatrix matsliceC;  //popB males
NumericMatrix matsliceD;  //popB females
for(int i=0;i<ages;i++){
matsliceC=NumericMatrix(N_init_a2_FIRST[i],nSNP);
matsliceD=NumericMatrix(N_init_a2_FIRST[i],nSNP);
for (int j=std::accumulate(N_init_a2_FIRST.begin(),N_init_a2_FIRST.begin()+i,0.0);
j<std::accumulate(N_init_a2_FIRST.begin(),N_init_a2_FIRST.begin()+i+1,0.0);j++){
matsliceC(j-std::accumulate(N_init_a2_FIRST.begin(),N_init_a2_FIRST.begin()+i,0.0),_)=allelesB(j,_);}
for (int j=std::accumulate(N_init_a2_FIRST.begin(),N_init_a2_FIRST.begin()+i,0.0)+pop2_size;
j<std::accumulate(N_init_a2_FIRST.begin(),N_init_a2_FIRST.begin()+i+1,0.0)+pop2_size;j++){
matsliceD(j-std::accumulate(N_init_a2_FIRST.begin(),N_init_a2_FIRST.begin()+i,0.0)-pop2_size,_)=allelesB(j,_);}

List temp6(2); temp6[0]=matsliceC;temp6[1]=matsliceD;
popB[ages-1-i]=temp6;
}

NumericVector popA_first21 = rev(N_init_a1_FIRST);
NumericVector popB_first21 = rev(N_init_a2_FIRST);
NumericVector newrecsA = allA.column(ages-1);
NumericVector newrecsB = allB.column(ages-1);
NumericVector dimspopA(Nyrs+ages);
NumericVector dimspopB(Nyrs+ages);
NumericVector pluscounterA_list(Nyrs,-2.0);
NumericVector pluscounterB_list(Nyrs,-2.0);
for(int i=0;i<ages;i++){
dimspopA[i]=popA_first21[i];
dimspopB[i]=popB_first21[i];
}
for(int i=ages;i<(Nyrs+ages+1);i++){
dimspopA[i]=newrecsA[i-ages+1];
dimspopB[i]=newrecsB[i-ages+1];
}

for(int k=0;k<Nyrs;k++){    //Nyrs

//set up spawning groups (for males and females, popA)***************************************************
Rcout<<"k at 2079"<<std::endl;
Rf_PrintValue(wrap(k));

spawn_malA_thisyr = NumericMatrix(2*NspawnersA[k],nSNP); //twice as many genotypes as spawners
spawn_femA_thisyr = NumericMatrix(2*NspawnersA[k],nSNP);
some_integerA=0;
spawnersA_rounded_star = allspawnA.row(k);  //from 21 yr olds to newborns
NumericVector N_init_a1_star = allA.row(k);  //also from 21 yr olds to newborns first row is N_init_a1_FIRST
tmpMatrixA_1 = NumericMatrix(0,nSNP);
tmpMatrixA_2 = NumericMatrix(0,nSNP);

//fix plus group problems
if ((spawnersA_rounded_star[0]>0)){ //only worry about them if there is a plus group
IntegerVector vec_plusA(k+1);
int plus_counterA = -1;
for(int j=0;j<k+1;j++){
if(j==0){
vec_plusA[j]=dimspopA(k-j);}  //vec_plus is a vector of the cumulative number of rows in popA from the current plus group to previous years, etc.
else{
vec_plusA[j]=vec_plusA[j-1]+dimspopA(k-j);}
if (vec_plusA[j]<2*spawnersA_rounded_star[0]) //plus_counterA tells you how far back you need to go
{plus_counterA = j;}
}
pluscounterA_list[k]=plus_counterA;
for (int j=k-(plus_counterA+1); j<k+1;j++){
tmpList_plusgroupA = popA[j];  //this is the actual plus group
tmpMatrixA_1plus=as<SEXP>(tmpList_plusgroupA[0]);
tmpMatrixA_2plus=as<SEXP>(tmpList_plusgroupA[1]);
for (int l=0;l<tmpMatrixA_1plus.nrow();l++){  
spawn_malA_thisyr(some_integerA,_)=tmpMatrixA_1plus(l,_);
spawn_femA_thisyr(some_integerA,_)=tmpMatrixA_2plus(l,_);
some_integerA = some_integerA+1;
if(some_integerA>=2*spawnersA_rounded_star[0]){j=k+1;l=tmpMatrixA_1plus.nrow();}
}
}

}//end of if spawnersA_rounded_star[0]>0  

//now your spawn_femA_thisyr has plus group in it.
//use some_integerA to fill in the rest of spawn_femA_thisyr
for(int i=1;i<ages;i++){                                   
if((spawnersA_rounded_star[i]>0)){
mal_spawnersA = IntegerVector(spawnersA_rounded_star[i]); //after random mating  **also make sure you get the right row
fem_spawnersA=IntegerVector(spawnersA_rounded_star[i]); //after random mating//may be 0 if spawnersA_rounded_star is 1 b/c of indexing (0 is 1)
spawner_listM = seq_len(N_init_a1_star[i]/2); //list of ones to choose from
spawner_listF = seq_len(N_init_a1_star[i]/2); //choosing from half the number of individuals because also using the subsequent row
tmpMatrixA_1 = NumericMatrix(0,nSNP);  
tmpMatrixA_2 = NumericMatrix(0,nSNP);

tmpListA = popA[i+k];  
tmpMatrixA_1=as<SEXP>(tmpListA[0]);
tmpMatrixA_2=as<SEXP>(tmpListA[1]);


//if (tmpMatrixA_1.nrow()>(2*N_init_a1_star[i])){

for (int j=0; j < spawnersA_rounded_star[i]; j++){ //random mating for each allele
pick_spawner_mal=as<double>(runif(1,-0.5,((N_init_a1_star[i]/2)-j-0.5)));//b/c choosing spawnersA_rounded_star[i] from all available and from half so can use subsequent row as well
pickM=int(pick_spawner_mal);
if((pick_spawner_mal+0.5)>=(int(pick_spawner_mal)+1)){  
pickM=int(pick_spawner_mal)+1;   
}
pick_spawner_fem=as<double>(runif(1,-0.5,((N_init_a1_star[i]/2)-j-0.5)));//b/c you are choosing spawnersA_rounded[i] from all those available
pickF=int(pick_spawner_fem);
if((pick_spawner_fem+0.5)>=(int(pick_spawner_fem)+1)){  
pickF=int(pick_spawner_fem)+1;   
}

mal_spawnersA[j]=spawner_listM[pickM]; //numeric(0) happens when you skip one because it is zero.
fem_spawnersA[j]=spawner_listF[pickF];  
spawn_malA_thisyr(some_integerA,_)=tmpMatrixA_1(2*spawner_listM[pickM]-2,_);  //the minus two is to get the indices to work out
spawn_malA_thisyr(some_integerA+1,_)=tmpMatrixA_1(2*spawner_listM[pickM]-1,_);//and the next one as well in the row
spawn_femA_thisyr(some_integerA,_)=tmpMatrixA_2(2*spawner_listF[pickF]-2,_);
spawn_femA_thisyr(some_integerA+1,_)=tmpMatrixA_2(2*spawner_listF[pickF]-1,_);
spawner_listM.erase(pickM);  //this erases the one at the location you specify
spawner_listF.erase(pickF);
some_integerA=some_integerA+2;

}  //end js for A
//}//end if loop

} //end if spawnersA_rounded_star[i]>0

}  //end of i loop (done with finding spawners)

//set up spawning groups (for males and females, popB)***************************************************
spawn_malB_thisyr = NumericMatrix(2*NspawnersB[k],nSNP); //twice as many genotypes as spawners
spawn_femB_thisyr = NumericMatrix(2*NspawnersB[k],nSNP);
some_integerB=0;
spawnersB_rounded_star = allspawnB.row(k);  //from 21 yr olds to newborns
NumericVector N_init_a2_star = allB.row(k);  //also from 21 yr olds to newborns first row is N_init_a1_FIRST
tmpMatrixB_1 = NumericMatrix(0,nSNP);
tmpMatrixB_2 = NumericMatrix(0,nSNP);
if (spawnersB_rounded_star[0]>0){ //only worry about them if there is a plus group
IntegerVector vec_plusB(k+1);
int plus_counterB = -1;
for(int j=0;j<k+1;j++){
if(j==0){
vec_plusB[j]=dimspopB(k-j);}  //vec_plus is a vector of the cumulative number of rows in popA from the current plus group to previous years, etc.
else{
vec_plusB[j]=vec_plusB[j-1]+dimspopB(k-j);}
if (vec_plusB[j]<2*spawnersB_rounded_star[0]) //plus_counterB tells you how far back you need to go
{plus_counterB = j;}
}
pluscounterB_list[k]=plus_counterB;
//problem is that by default it is 0 so it does not help 
for (int j=k-(plus_counterB+1); j<k+1;j++){
tmpList_plusgroupB = popB[j];  //this is the actual plus group
tmpMatrixB_1plus=as<SEXP>(tmpList_plusgroupB[0]);
tmpMatrixB_2plus=as<SEXP>(tmpList_plusgroupB[1]);
for (int l=0;l<tmpMatrixB_1plus.nrow();l++){  
spawn_malB_thisyr(some_integerB,_)=tmpMatrixB_1plus(l,_);
spawn_femB_thisyr(some_integerB,_)=tmpMatrixB_2plus(l,_);
some_integerB = some_integerB+1;
if(some_integerB>=2*spawnersB_rounded_star[0]){j=k+1;l=tmpMatrixB_1plus.nrow();}  //>=because of C++ indexing
}
}

}//end of if spawnersB_rounded_star[0]>0  

//now your spawn_femB_thisyr has plus group in it.

for(int i=1;i<ages;i++){                                  
if((spawnersB_rounded_star[i]>0)){
mal_spawnersB = IntegerVector(spawnersB_rounded_star[i]); //after random mating  **also make sure you get the right row
fem_spawnersB = IntegerVector(spawnersB_rounded_star[i]); //after random mating
spawner_listM = seq_len(N_init_a2_star[i]/2); //list of ones to choose from
spawner_listF = seq_len(N_init_a2_star[i]/2); //choosing from half the number of individuals because also using the subsequent row
tmpMatrixB_1 = NumericMatrix(0,nSNP);  
tmpMatrixB_2 = NumericMatrix(0,nSNP); 

tmpListB = popB[i+k];  
tmpMatrixB_1=as<SEXP>(tmpListB[0]);
tmpMatrixB_2=as<SEXP>(tmpListB[1]);

for (int j=0; j < spawnersB_rounded_star[i]; j++){ //random mating for each allele
pick_spawner_mal=as<double>(runif(1,-0.5,((N_init_a2_star[i]/2)-j-0.5)));//b/c choosing spawnersB_rounded_star[i] from all available and from half so can use subsequent row as well
pickM=int(pick_spawner_mal);
if((pick_spawner_mal+0.5)>=(int(pick_spawner_mal)+1)){  
pickM=int(pick_spawner_mal)+1;   
}

pick_spawner_fem=as<double>(runif(1,-0.5,((N_init_a2_star[i]/2)-j-0.5)));//b/c you are choosing spawnersA_rounded[i] from all those available
pickF=int(pick_spawner_fem);
if((pick_spawner_fem+0.5)>=(int(pick_spawner_fem)+1)){  
pickF=int(pick_spawner_fem)+1;   
}

mal_spawnersB[j]=spawner_listM[pickM]; //numeric(0) happens when you skip one because it is zero.
fem_spawnersB[j]=spawner_listF[pickF];  

spawn_malB_thisyr(some_integerB,_)=tmpMatrixB_1(2*spawner_listM[pickM]-2,_);  //the minus two is to get the indices to work out
spawn_malB_thisyr(some_integerB+1,_)=tmpMatrixB_1(2*spawner_listM[pickM]-1,_);//and the next one as well in the row

spawn_femB_thisyr(some_integerB,_)=tmpMatrixB_2(2*spawner_listF[pickF]-2,_);
spawn_femB_thisyr(some_integerB+1,_)=tmpMatrixB_2(2*spawner_listF[pickF]-1,_);

spawner_listM.erase(pickM);  //this erases the one at the location you specify
spawner_listF.erase(pickF);

some_integerB=some_integerB+2;

}  //end js for B
} //end if spawnersB_rounded_star[i]>0

}  //end of i loop (done with finding spawners)


spawn_femA[k]=spawn_femA_thisyr;
spawn_malA[k]=spawn_malA_thisyr;
spawn_femB[k]=spawn_femB_thisyr;
spawn_malB[k]=spawn_malB_thisyr;

//now establish genotypes for the new recruits
RecA = NumericMatrix (2*allA(k+1,ages-1),nSNP);
RecB = NumericMatrix (2*allB(k+1,ages-1),nSNP);
//recruits for popA

for(int j=0;j<allA(k+1,ages-1);j++){
rand_doub=as<double>(runif(1,-0.5,(spawn_femA_thisyr.nrow()/2)-1));
samp_mom=int(rand_doub);
if((rand_doub+0.5)>=(int(rand_doub)+1)){  
samp_mom=int(rand_doub)+1;   
}
rand_doub=as<double>(runif(1,-0.5,(spawn_malA_thisyr.nrow()/2)-1));
samp_dad=int(rand_doub);
if((rand_doub+0.5)>=(int(rand_doub)+1)){	
samp_dad=int(rand_doub)+1;   
}

for(int i=0;i<nSNP;i++){
rand_doub=as<double>(runif(1,-0.5,1.5));
rand_allele=int(rand_doub);
if((rand_doub+0.5)>=(int(rand_doub)+1)){  
rand_allele=int(rand_doub)+1;   
}
RecA(2*j,i)=spawn_malA_thisyr(samp_dad*2+rand_allele,i);
RecA(2*j+1,i)=spawn_femA_thisyr(samp_mom*2+rand_allele,i);

}  //end of all 10 microsatellites

}  //end of j RecAs recruits popA

//recruits for popB
for(int j=0;j<allB(k+1,ages-1);j++){
rand_doub=as<double>(runif(1,-0.5,(spawn_femB_thisyr.nrow()/2)-1));
samp_mom=int(rand_doub);
if((rand_doub+0.5)>=(int(rand_doub)+1)){  
samp_mom=int(rand_doub)+1;   
}
rand_doub=as<double>(runif(1,-0.5,(spawn_malB_thisyr.nrow()/2)-1));
samp_dad=int(rand_doub);
if((rand_doub+0.5)>=(int(rand_doub)+1)){	
samp_dad=int(rand_doub)+1;   
}
for(int i=0;i<nSNP;i++){
rand_doub=as<double>(runif(1,-0.5,1.5));
rand_allele=int(rand_doub);
if((rand_doub+0.5)>=(int(rand_doub)+1)){  
rand_allele=int(rand_doub)+1;   
}
RecB(2*j,i)=spawn_malB_thisyr(samp_dad*2+rand_allele,i);
RecB(2*j+1,i)=spawn_femB_thisyr(samp_mom*2+rand_allele,i);
}  //end of all 10 microsatellites
}  //end of j RecBs recruits popB
//here take the matrices RecA and RecB which are random offspring
//and select a subset of the individuals and make more of them
//this represents a subset of parents being successful at mating
//there might be better ways, like reducing number of spawners

recA = clone(RecA);
recB = clone(RecB);

//mutate new recruits
double prob;
double dir;
for (int j=0;j<nSNP;j++){
for (int i=0; i < recA.nrow(); i++){  //recA mutation for each allele
prob=as<double>(runif(1,0,1));  
dir=as<double>(runif(1,0,1));		
if(prob<0.01*MUS[j]&dir<0.5){	
recA(i,j)=recA(i,j)+1;
}else if (prob<MUS[j]&dir>0.5){
recA(i,j)=recA(i,j)-1;
}			
}
for (int i=0; i < recB.nrow(); i++){  //recB mutation for each allele
prob=as<double>(runif(1,0,1));  
dir=as<double>(runif(1,0,1));		
if(prob<0.01*MUS[j]&dir<0.5){	
recB(i,j)=recB(i,j)+1;
}else if (prob<MUS[j]&dir>0.5){
recB(i,j)=recB(i,j)-1;
}			
}
}

//MIGRATION


if ((k>1)&(k%1==0)){
Rcout<<"second migration"<<std::endl;
nmig=nmig1;


NumericVector BtoAlist(nmig);
NumericVector AfromBlist(nmig);
int whichone;
if (recA.nrow()<=recB.nrow()){whichone=1;}
if (recB.nrow()<recA.nrow()){whichone=2;}
if (whichone==1){//if there are fewer recruits in recA
if ((2*nmig)>recA.nrow()|(2*nmig)>recB.nrow()){nmig=recA.nrow()/2;}
}
if (whichone==2){//if there are fewer recruits in recB
if ((2*nmig)>recA.nrow()|(2*nmig)>recB.nrow()){nmig=recB.nrow()/2;}
}

Rcout<<"nmig,recA,recB,whichone"<<std::endl;
Rf_PrintValue(wrap(nmig));
Rf_PrintValue(wrap(recA.nrow()));
Rf_PrintValue(wrap(recB.nrow()));
Rf_PrintValue(wrap(whichone));



nmig_vec[k]=nmig;


NumericMatrix clonerecB=clone(recB);
NumericMatrix clonerecA=clone(recA);

BtoAlist=SeqLen(clonerecB.nrow()/2)-1;
AfromBlist=SeqLen(clonerecA.nrow()/2)-1;

for (int i=0;i<nmig;i++){
double probfromB=as<double>(runif(1,0,(clonerecB.nrow()/2)-1-i)); //works, I checked it
double probtoA=as<double>(runif(1,0,(clonerecA.nrow()/2)-1-i));

int pick_fromB=int(probfromB);
if((probfromB+0.5)>=(int(probfromB)+1)){
pick_fromB=int(probfromB)+1;
}

int pick_toA=int(probtoA);
if((probtoA+0.5)>=(int(probtoA)+1)){
pick_toA=int(probtoA)+1;
}

 NumericMatrix::Row moverB2Arow1=clonerecB(2*(BtoAlist[pick_fromB]),_);  //minus 1 is for indexing
 NumericMatrix::Row moverB2Arow2=clonerecB((2*(BtoAlist[pick_fromB]))+1,_);

  recA(2*(AfromBlist[pick_toA]),_)=moverB2Arow1;
  recA((2*(AfromBlist[pick_toA]))+1,_)=moverB2Arow2;

BtoAlist.erase(pick_fromB);
AfromBlist.erase(pick_toA);
}//end of migration

//for (int i=0;i<2*nmig;i++){

//  NumericMatrix::Row moverB2Arow=clonerecB(i,_);  //minus 1 is for indexing
//  recA(i,_)=moverB2Arow;
//}


}//end of if k >1


//save new recruits to popA and popB
NumericMatrix recA_1(allA(k+1,ages-1),nSNP);
NumericMatrix recA_2(allA(k+1,ages-1),nSNP);
for(int i=0;i<allA(k+1,ages-1);i++){
recA_1(i,_)=recA(i,_);
recA_2(i,_)=recA(i+allA(k+1,ages-1),_);
}

tmpList_recsA[0] = recA_1;
tmpList_recsA[1] = recA_2;
NumericMatrix recB_1(allB(k+1,ages-1),nSNP);
NumericMatrix recB_2(allB(k+1,ages-1),nSNP);
for(int i=0;i<allB(k+1,ages-1);i++){
recB_1(i,_)=recB(i,_);
recB_2(i,_)=recB(i+allB(k+1,ages-1),_);
}
tmpList_recsB[0] = recB_1;
tmpList_recsB[1] = recB_2;
popA[k+ages]=clone(tmpList_recsA);
popB[k+ages]=clone(tmpList_recsB);

}// end k loop

//get genetic diversity - put all genotypes in GenotypesA1, GenotypesA2, GenotypesB1, GenotypesB2
NumericVector cumpopsA(Nyrs+ages);//cumpops holds the cumulative number (for counting purposes)
NumericVector cumpopsB(Nyrs+ages);
cumpopsA[0]=0;
cumpopsB[0]=0;
for(int i=1;i<(Nyrs+ages);i++){
cumpopsA[i]=dimspopA[i-1]+cumpopsA[i-1];
cumpopsB[i]=dimspopB[i-1]+cumpopsB[i-1];
}
NumericMatrix GenotypesA1(cumpopsA[Nyrs+ages-1],nSNP);//matrix to hold everything  (we want this but come back to it)
NumericMatrix GenotypesA2(cumpopsA[Nyrs+ages-1],nSNP);
NumericMatrix GenotypesB1(cumpopsB[Nyrs+ages-1],nSNP);
NumericMatrix GenotypesB2(cumpopsB[Nyrs+ages-1],nSNP);
NumericMatrix tmpMatA_1;
NumericMatrix tmpMatA_2;
NumericMatrix tmpMatB_1;
NumericMatrix tmpMatB_2;
List tmpLA(2);
List tmpLB(2);
some_integerA=0;
some_integerB=0;
Rf_PrintValue(wrap(1714));
//pluscounterB_list is Nyrs long (110)

NumericVector poprowsA(Nyrs+ages-1);
for (int i=0;i<(Nyrs+ages-1);i++){
List hmmA(2);
hmmA=popA[i];
NumericMatrix hmm2A=hmmA[1];
poprowsA[i]=hmm2A.nrow();}
GenotypesA1=NumericMatrix(std::accumulate(poprowsA.begin(),poprowsA.end(),0.0),nSNP);
GenotypesA2=NumericMatrix(std::accumulate(poprowsA.begin(),poprowsA.end(),0.0),nSNP);

for(int k=0; k<(Nyrs+ages-1);k++){
tmpLA = popA[k];
tmpMatA_1 = as<SEXP>(tmpLA[0]);
tmpMatA_2 = as<SEXP>(tmpLA[1]);
for(int i=0; i<tmpMatA_1.nrow();i++){
GenotypesA1(some_integerA,_)=tmpMatA_1(i,_);
GenotypesA2(some_integerA,_)=tmpMatA_2(i,_);
some_integerA=some_integerA+1;
}
}

//genotypesB1<sum of all popB or <tmpMatB_
NumericVector poprowsB(Nyrs+ages-1);
for (int i=0;i<(Nyrs+ages-1);i++){
List hmm(2);
hmm=popB[i];
NumericMatrix hmm2=hmm[1];
poprowsB[i]=hmm2.nrow();}
GenotypesB1=NumericMatrix(std::accumulate(poprowsB.begin(),poprowsB.end(),0.0),nSNP);
GenotypesB2=NumericMatrix(std::accumulate(poprowsB.begin(),poprowsB.end(),0.0),nSNP);

//if (genotypesB1.nrow()<
for(int k=0; k<(Nyrs+ages-1);k++){
tmpLB = popB[k];
tmpMatB_1 = as<SEXP>(tmpLB[0]);
tmpMatB_2 = as<SEXP>(tmpLB[1]);
for(int i=0; i<tmpMatB_1.nrow();i++){
GenotypesB1(some_integerB,_)=tmpMatB_1(i,_);
GenotypesB2(some_integerB,_)=tmpMatB_2(i,_);
some_integerB=some_integerB+1;
}
Rcout<<"k in this loop goes to 117 but 117 is empty"<<std::endl;
Rf_PrintValue(wrap(k));
}


//get allelic richness for popA - the initial year
all_richA=NumericVector(nSNP);  //set up a vector to hold number of alleles for each microsat
int rich_nA = arichAmat.nrow();
for (int rich = 0; rich<nSNP;rich++){
NumericMatrix::Column inputA = arichAmat(_,rich); 
NumericVector xA = clone<NumericVector>(inputA);

int nallsA=0;
typedef std::map<double,int> imap ;
typedef imap::value_type pair ;
imap index ;
int n = xA.size() ;
double current, previous = xA[0] ;
index.insert( pair( previous, 0 ) );
imap::iterator it = index.begin() ;
for( int i=1; i<n; i++){
current = xA[i] ;
if( current == previous ){
xA[i] = current + ( ++(it->second) / 100.0 ) ;
} else {
it = index.find(current) ;
if( it == index.end() ){
it = index.insert(
current > previous ? it : index.begin(),
pair( current, 0 )
) ;
} else {
xA[i] = current + ( ++(it->second) / 100.0 ) ;
}
previous = current ;
}

if (xA[i]-inputA[i]==0){nallsA=nallsA+1;}
}
all_richA[rich] = nallsA+1;
} //end of rich A loop (over 10 microsats)
//end of rch code popA
richnessA(0,_)=all_richA;
mean_arichA[0]=std::accumulate(all_richA.begin(),all_richA.end(),0.0)/nSNP;

//get allelic richness for popB - the initial year
all_richB=NumericVector(nSNP);  //set up a vector to hold number of alleles for each microsat
int rich_nB = arichBmat.nrow();
for (int rich = 0; rich<nSNP;rich++){
NumericMatrix::Column inputB = arichBmat(_,rich); 
NumericVector xB = clone<NumericVector>(inputB);
int nallsB=0;
typedef std::map<double,int> imap ;
typedef imap::value_type pair ;
imap index ;
int n = xB.size() ;
double current, previous = xB[0] ;
index.insert( pair( previous, 0 ) );
imap::iterator it = index.begin() ;
for( int i=1; i<n; i++){
current = xB[i] ;
if( current == previous ){
xB[i] = current + ( ++(it->second) / 100.0 ) ;
} else {
it = index.find(current) ;
if( it == index.end() ){
it = index.insert(
current > previous ? it : index.begin(),
pair( current, 0 )
) ;
} else {
xB[i] = current + ( ++(it->second) / 100.0 ) ;
}
previous = current ;
}

if (xB[i]-inputB[i]==0){nallsB=nallsB+1;}
}
all_richB[rich] = nallsB+1;
} //end of rich B loop (over 10 microsats)
//end of rch code popB
richnessB(0,_)=all_richB;  //each row is allelic richness at each microsat in that year
mean_arichB[0]=std::accumulate(all_richB.begin(),all_richB.end(),0.0)/nSNP;

//get allelic richness in other years popA
for (int k = 1; k<Nyrs; k++){  //Nyrs
some_integerA = 0;
arichAmat = NumericMatrix(2*NtotalA[k],nSNP);
N_init_a1_star = allA.row(k);

int j_a; int j_b;

for (int i=0;i<ages;i++){
if ((cumpopsA[k+i]+N_init_a1_star[i])<=cumpopsA[Nyrs+ages-1]){
j_a = cumpopsA[k+i]; j_b = cumpopsA[k+i]+N_init_a1_star[i];}else
{j_a = cumpopsA[Nyrs+ages-1]-N_init_a1_star[i]; j_b = cumpopsA[Nyrs+ages-1];}

for(int j=j_a;j<j_b;j++){
arichAmat(some_integerA,_)=GenotypesA1(j,_);
arichAmat(some_integerA+NtotalA[k],_)=GenotypesA2(j,_);
some_integerA=some_integerA+1;
}
}

//get allelic richness for popA (all years except the first one)
all_richA=NumericVector(nSNP);  //set up a vector to hold number of alleles for each microsat
int rich_nA = arichAmat.nrow();
for (int rich = 0; rich<nSNP;rich++){
NumericMatrix::Column inputA = arichAmat(_,rich); 
NumericVector xA = clone<NumericVector>(inputA);  //already defined xA

int nallsA=0;
typedef std::map<double,int> imap ;
typedef imap::value_type pair ;
imap index ;
int n = xA.size() ;
double current, previous = xA[0] ;
index.insert( pair( previous, 0 ) );
imap::iterator it = index.begin() ;
for( int i=1; i<n; i++){
current = xA[i] ;
if( current == previous ){
xA[i] = current + ( ++(it->second) / 100.0 ) ;
} else {
it = index.find(current) ;
if( it == index.end() ){
it = index.insert(
current > previous ? it : index.begin(),
pair( current, 0 )
) ;
} else {
xA[i] = current + ( ++(it->second) / 100.0 ) ;
}
previous = current ;
}

if (xA[i]-inputA[i]==0)
{nallsA=nallsA+1;
}
}
//*****************************
all_richA[rich] = nallsA+1;
} //end of rich A loop (over 10 microsats)
//end of rch code popA
richnessA(k,_)=all_richA;
mean_arichA[k]=std::accumulate(all_richA.begin(),all_richA.end(),0.0)/nSNP;
//		}//end of k= 1 to Nyrs popA  
//get allelic richness in other years popB
//		for (int k = 1; k<Nyrs; k++){  //Nyrs
some_integerB = 0;
arichBmat = NumericMatrix(2*NtotalB[k],nSNP);
N_init_a2_star = allB.row(k);
//int j_a; int j_b; 
for (int i=0;i<ages;i++){

if ((cumpopsB[k+i]+N_init_a2_star[i])<=cumpopsB[Nyrs+ages-1]){
j_a = cumpopsB[k+i]; j_b = cumpopsB[k+i]+N_init_a2_star[i];}else
{j_a = cumpopsB[Nyrs+ages-1]-N_init_a2_star[i]; j_b = cumpopsB[Nyrs+ages-1];}

for(int j=j_a;j<j_b;j++){
arichBmat(some_integerB,_)=GenotypesB1(j,_);
arichBmat(some_integerB+NtotalB[k],_)=GenotypesB2(j,_);
some_integerB=some_integerB+1;
}
}
Rf_PrintValue(wrap(1658));
//get allelic richness for popB
all_richB=NumericVector(nSNP);  //set up a vector to hold number of alleles for each microsat
int rich_nB = arichBmat.nrow();
Rf_PrintValue(wrap(1665));
for (int rich = 0; rich<nSNP;rich++){
NumericMatrix::Column inputB = arichBmat(_,rich); 
NumericVector xB = clone<NumericVector>(inputB);

int nallsB=0;
typedef std::map<double,int> imap ;
typedef imap::value_type pair ;
imap index ;
int n = xB.size() ;
double current, previous = xB[0] ;
index.insert( pair( previous, 0 ) );
imap::iterator it = index.begin() ;
for( int i=1; i<n; i++){
current = xB[i] ;
if( current == previous ){
xB[i] = current + ( ++(it->second) / 100.0 ) ;
} else {
it = index.find(current) ;
if( it == index.end() ){
it = index.insert(
current > previous ? it : index.begin(),
pair( current, 0 )
) ;
} else {
xB[i] = current + ( ++(it->second) / 100.0 ) ;
}
previous = current ;
}

if (xB[i]-inputB[i]==0)
{nallsB=nallsB+1;
}
}

all_richB[rich] = nallsB+1;
} //end of rich B loop (over 10 microsats)

//end of rch code popB
richnessB(k,_)=all_richB;  //each row is allelic richness at each microsat in that year

mean_arichB[k]=std::accumulate(all_richB.begin(),all_richB.end(),0.0)/nSNP;


//for (int j=0;j<nSNP;j++){
//for (int i=0;i<(arichAmat.nrow()/2);i++){
//if (arichAmat((2*i),j)==arichAmat(((2*i)+1),j)){hetsA[j]=hetsA[j]+1;}
//}}

//for (int j=0;j<nSNP;j++){
//for (int i=0;i<(arichBmat.nrow()/2);i++){
//if (arichBmat((2*i),j)==arichBmat(((2*i)+1),j)){hetsB[j]=hetsB[j]+1;}
//}}

//***********************************************************************
//here get new FST based on this years arichAmat and arichBmat
//arichAmat and arichBmat get generated each of k years so they are already available.
//Fst (for the entire population)  if nsamp is bigger than pop size, will give warning
//create a matrix with unique alleles for each microsatellite in each column

if ((k%20==0)&(NtotalA[k]>0)&(NtotalB[k]>0)&(mgmt_counter==0)){   //do this for all years for combined mgmt
//for all k
NumericMatrix tallyBOTHmat(2*NtotalB[k]+2*NtotalA[k],nSNP);
int some_integerBOTH = 0;
NumericMatrix arichBOTHmat(2*NtotalB[k]+2*NtotalA[k],nSNP);//all alleles both pops year 1

for (int i=0;i<2*NtotalA[k];i++){
arichBOTHmat(i,_)=arichAmat(i,_);
}
for (int i=0;i<2*NtotalB[k];i++){
arichBOTHmat(i+2*NtotalA[k],_)=arichBmat(i,_);
}
//get unique alleles for both populations
for (int rich = 0; rich<nSNP;rich++){
NumericMatrix::Column inputBOTH = arichBOTHmat(_,rich); 
NumericVector xBOTH = clone<NumericVector>(inputBOTH);
NumericVector tally_allelesBOTH(2*NtotalB[k]+2*NtotalA[k]); //this is a vector that holds the unique alleles
int tallycounter=1;
//*******
int nallsBOTH=0;
typedef std::map<double,int> imap ;
typedef imap::value_type pair ;
imap index ;
int n = xBOTH.size() ;
double current, previous = xBOTH[0] ;
index.insert( pair( previous, 0 ) );
imap::iterator it = index.begin() ;
for( int i=1; i<n; i++){
current = xBOTH[i] ;
if( current == previous ){
xBOTH[i] = current + ( ++(it->second) / 100.0 ) ;
} else {
it = index.find(current) ;
if( it == index.end() ){
it = index.insert(
current > previous ? it : index.begin(),
pair( current, 0 )
) ;
} else {
xBOTH[i] = current + ( ++(it->second) / 100.0 ) ;
}
previous = current ;
}

if (xBOTH[i]-inputBOTH[i]==0)
{nallsBOTH=nallsBOTH+1;
tally_allelesBOTH[tallycounter]=inputBOTH[i];
tallycounter=tallycounter+1;
}}
tally_allelesBOTH[0]=inputBOTH[0];
tallyBOTHmat(_,rich)=tally_allelesBOTH;
} //end of rich BOTH loop (over 10 microsats)

List tallist(nSNP); //list of 10 that will hold allele frequencies for each of 10 microsatellites.
//First column is the actual allele size, second is freq from popA, third is freq popB
NumericVector count(nSNP);//list of number of unique alleles at each locus
//now set up a table to get allele frequencies for popA and popB
for (int rich=0;rich<nSNP;rich++){
NumericMatrix::Column inputA=arichAmat(_,rich);
NumericMatrix::Column inputB=arichBmat(_,rich);
NumericMatrix::Column tal=tallyBOTHmat(_,rich); 
for(int i=0;i<tal.size();i++){  //figure out how many unique alleles are at this microsat (from tal)
if(tal[i]!=0){count[rich]=count[rich]+1;}
}

NumericVector Counter=clone(count);
int counter=Counter[rich];
NumericVector taltrunc(counter);  //will hold all unique alleles from this microsat (both pops)
NumericMatrix tallyAB(counter,3);  //matrix that has alleles, freq at popA, freq at popB
NumericVector howmanyA(counter);  //number of alleles for this population at each allele
NumericVector howmanyB(counter);

for(int i=0;i<counter;i++){  //counter is the number of unique alleles at a locus
taltrunc[i]=tal[i];
int counterA=0;  //a counter for number of unique alleles at each locus
int counterB=0;

for (int j=0;j<2*NtotalA[k];j++){
if (inputA[j]==taltrunc[i])//go through all alleles to see how many match this unique one
{counterA=counterA+1;}
}
howmanyA[i]=counterA;
for (int j=0;j<2*NtotalB[k];j++){
if (inputB[j]==taltrunc[i])
{counterB=counterB+1;}
}
howmanyB[i]=counterB;

} //end of counter
tallyAB(_,0)=taltrunc;
tallyAB(_,1)=howmanyA/(2*NtotalA[k]);
tallyAB(_,2)=howmanyB/(2*NtotalB[k]);
tallist[rich]=tallyAB;

}//end of rich create talmat, which has unique alleles first column then freqs for popA and B
NumericMatrix talmat(std::accumulate(count.begin(),count.end(),0.0),5);

for (int i=0;i<nSNP;i++){
int talcount=0;
NumericMatrix taltmp(count[i],3);
taltmp=as<SEXP>(tallist[i]);
for (int j=std::accumulate(count.begin(),count.begin()+i,0.0);j<std::accumulate(count.begin(),count.begin()+i+1,0.0);j++)
{
talmat(j,_)=taltmp(talcount,_);
talcount=talcount+1;
}
}


double hetsA;
double hetsB;
for (int i=0;i<nSNP;i++){

for (int j=std::accumulate(count.begin(),count.begin()+i,0.0);j<std::accumulate(count.begin(),count.begin()+i+1,0.0);j++)
{
hetsA=0;hetsB=0;
for (int k=0;k<100;k++){
if ((arichAmat((2*k),i)!=arichAmat(((2*k)+1),i))&((arichAmat((2*k),i)==talmat(j,0))|(arichAmat(((2*k)+1),i)==talmat(j,0)))){hetsA=hetsA+1;}
if ((arichBmat((2*k),i)!=arichBmat(((2*k)+1),i))&((arichBmat((2*k),i)==talmat(j,0))|(arichBmat(((2*k)+1),i)==talmat(j,0)))){hetsB=hetsB+1;}
}
talmat(j,3)=hetsA/pop1_size;
talmat(j,4)=hetsB/pop2_size;
}}   

talmat2=clone(talmat);


//talmat is the raw material to find Fst!!!
//GET FST this function finds Fst between 2 populations, multiple loci and multiple alleles.
n_sampA=NtotalA[k];
n_sampB=NtotalB[k];
double n_bar=0.5*(n_sampA+n_sampB); //sample size - can change this.
double r=2;

NumericVector p_bar(talmat2.nrow());
NumericMatrix::Column sampmatA =talmat2(_,1);
NumericMatrix::Column sampmatB =talmat2(_,2);
NumericVector s2(talmat2.nrow());
NumericVector h_bar(talmat2.nrow());
NumericVector ones(talmat2.nrow());
ones=rep(1,talmat2.nrow());
p_bar=(n_sampA*sampmatA+n_sampB*sampmatB)/(2*n_bar);  //each entry is average sample frequency of an allele
s2=((n_sampB*(pow((sampmatB-p_bar),2)))+(n_sampA*pow((sampmatA-p_bar),2)))/n_bar;//sample variance of allele freqs over pops (for each allele)

for (int i=0;i<talmat2.nrow();i++){
h_bar[i]=(talmat2(i,3)+talmat2(i,4))/2;
}

double nc=((r*n_bar)-((pow(n_sampA,2)/(r*n_bar))+(pow(n_sampB,2)/(r*n_bar))))/(r-1);  //same as n_bar

double C2=r*(1-(nc/n_bar));
Rcout<<"C2"<<std::endl;
Rf_PrintValue(wrap(C2));

NumericVector a=(n_bar/nc)*(s2-((1/(n_bar-1))*((p_bar*(1-p_bar))-(((r-1)/r)*s2)-(0.25*h_bar))));
NumericVector dvec=((2*n_bar)/((2*n_bar)-1))*((p_bar*(1-p_bar))-(((r-1)/r)*s2));
NumericVector b=(n_bar/(n_bar-1))*((p_bar*(1-p_bar))-(((r-1)/r)*s2)-(h_bar*(2*n_bar-1)/(4*n_bar)));
NumericVector c=h_bar/2;
NumericVector aplusdvec=a+b+c;//I added b here
NumericVector aplusdvec2=a+dvec;
double fst=std::accumulate(a.begin(),a.end(),0.0)/(std::accumulate(aplusdvec.begin(),aplusdvec.end(),0.0)); //the better one
double fst2=std::accumulate(a.begin(),a.end(),0.0)/(std::accumulate(aplusdvec2.begin(),aplusdvec2.end(),0.0)); //the better one

ps=p_bar*(1-p_bar);
NumericVector theta=(s2-((1/((2*n_bar)-1))*((p_bar*(1-p_bar))-s2)))/(((1-((2*n_bar*C2)/(((2*n_bar)-1)*r)))*(p_bar*(1-p_bar)))+(1+((2*n_bar*C2)/((2*n_bar-1)*r)))*(s2/r));

NumericVector theta1=(s2-((1/(n_bar-1))*((p_bar*(1-p_bar))-s2-(h_bar/4))))/((((1-((n_bar*C2)/((n_bar-1)*r)))*(p_bar*(1-p_bar)))+(((1+n_bar*C2)/((n_bar-1)*r))*(s2/r))+((C2/(r*(n_bar-1)))*(h_bar/4))));

out_fst60[0]=fst;
fst_vec[ngen+k]=fst;
fst_vec2[ngen+k]=fst2;
Wrights_fst_vec[ngen+k]=std::accumulate(theta.begin(),theta.end(),0.0)/talmat2.nrow();
Wrights_fst_vec1[ngen+k]=std::accumulate(theta1.begin(),theta1.end(),0.0)/talmat2.nrow();
Wrights_simple[ngen+k]=std::accumulate(s2.begin(),s2.end(),0.0)/std::accumulate(ps.begin(),ps.end(),0.0);

}//end of if statement
Rcout<<"k is"<<std::endl;

Rf_PrintValue(wrap(k));
//****GET FST FOR JUST A RANDOM SAMPLE OF 100
//if((k==10|k==20)){
if((k==620|k==640|k==680)&mgmt_counter==0){  //ends at 3106
n_sampA=nsamp;//just 100 generally 
n_sampB=nsamp;//and 100
IntegerVector order(2*n_sampA+2*n_sampB);
int nsig=100;
NumericVector fst_subvec(nsig);
IntegerVector ORDER(4*nsamp);
Rf_PrintValue(wrap(2186));


//randomly select 2*diffsamp without replacement from populationA
IntegerVector randsampA = seq_len(2*NtotalA[k]);
IntegerVector whichsampA(2*n_sampA);
NumericMatrix arichAmat_sub(2*n_sampA,nSNP);
NumericMatrix arichBmat_sub(2*n_sampB,nSNP);
int probA; int pickA;
Rf_PrintValue(wrap(2193));
if((NtotalA[k]>=n_sampA)&(NtotalB[k]>=n_sampB)){

for (int i=0; i < n_sampA; i++){  
probA=as<double>(runif(1,-0.5,(NtotalA[k]-0.5)));//sample without replacement
pickA=int(probA);
if((probA+0.5)>=(int(probA)+1)){  
pickA=int(probA)+1;
}
whichsampA[i]=randsampA[pickA]-1;//because of C++ special indexing
//randsampA.erase(pickA);
}

for (int i=0;i<n_sampA;i++){
arichAmat_sub(i,_)=arichAmat(whichsampA[i],_); //arichAmat is all alleles males and females from population that year
arichAmat_sub(i+n_sampA,_)=arichAmat(whichsampA[i]+1,_);
}  //arichAmat_sub is a random sample of 2*n_sampA
arichBmat_sub=NumericMatrix(2*n_sampB,nSNP);

//randomly select 2*diffsamp without replacement from populationB
//if (ninecomp[nine]==2){
IntegerVector randsampB = seq_len(2*NtotalB[k]);
IntegerVector whichsampB(2*n_sampB);
int probB;
int pickB;
for (int i=0; i < n_sampB; i++){
probB=as<double>(runif(1,-0.5,(NtotalB[k]-0.5)));//each time there is one fewer to choose from
pickB=int(probB);
if((probB+0.5)>=(int(probB)+1)){	
pickB=int(probB)+1;
}
whichsampB[i]=randsampB[pickB]-1;
//randsampB.erase(pickB);
}
for (int i=0;i<n_sampB;i++){
arichBmat_sub(i,_)=arichBmat(whichsampB[i],_);
arichBmat_sub(i+n_sampB,_)=arichBmat(whichsampB[i]+1,_);
}

//arichAmat_sub 

//} //end of if ninecomp[nine]==2;

if (k==620){arichAmat_sub620=arichAmat_sub;
arichBmat_sub620=arichBmat_sub;}
if (k==640){arichAmat_sub640=arichAmat_sub;
arichBmat_sub640=arichBmat_sub;}
if (k==680){arichAmat_sub680=arichAmat_sub;
arichBmat_sub680=arichBmat_sub;}



//by here you should have arichAmat_sub and arichBmat_sub create a matrix with unique alleles for each microsatellite in each column
NumericMatrix tallyBOTHmat_sub(2*n_sampA+2*n_sampB,nSNP); //this has everything subsampled, just 200+200 total
NumericMatrix arichBOTHmat_sub(2*n_sampA+2*n_sampB,nSNP); //all alleles both pops year 1
NumericMatrix arichBOTHmat_subrand(2*n_sampA+2*n_sampB,nSNP);
int some_integerBOTH_sub = 0;
for (int i=0;i<2*n_sampA;i++){    //holds 200 genotypes from each population
arichBOTHmat_sub(i,_)=arichAmat_sub(i,_);}
for (int i=0;i<2*n_sampB;i++){
arichBOTHmat_sub(i+2*n_sampA,_)=arichBmat_sub(i,_);
}
NumericVector fstsig(1); //this is going to count how many times the initial sample is larger than permuted ones
fstsig[0]=0;
for (int sigfst=0;sigfst<nsig;sigfst++){//want 100 later
ORDER=seq_len(2*n_sampA+2*n_sampB)-1;
order=rep(0,2*n_sampA+2*n_sampB);
if (sigfst==0)
{order=seq_len(2*n_sampA+2*n_sampB)-1;}
else{
int probZ;
int pickZ;
for (int i=0; i < 2*n_sampA+2*n_sampB; i++){  //random mating for each allele 4*nsamp
probZ=as<double>(runif(1,-0.5,(2*n_sampA+2*n_sampB-0.5-i)));
pickZ=int(probZ);
if((probZ+0.5)>=(int(probZ)+1)){	
pickZ=int(probZ)+1;
}

order[i]=ORDER[pickZ];
ORDER.erase(pickZ);
}
}//end of is this the first iteration of sigfst or not

for(int i=0;i<2*n_sampA+2*n_sampB;i++){  //rearrange arichBOTHmat_sub for significance calculations
arichBOTHmat_subrand(i,_)=arichBOTHmat_sub(order[i],_);
}
//divide up arichBOTHmat_subrand into 2 populations again for significance calculations first time these will be truly from pops A and B but later just jumbles of both
NumericMatrix arichAmat_subrand(2*n_sampA,nSNP);
NumericMatrix arichBmat_subrand(2*n_sampB,nSNP);
for (int i=0;i<2*n_sampA;i++){
arichAmat_subrand(i,_)=arichBOTHmat_subrand(i,_);}
for (int i=0;i<2*n_sampB;i++){
arichBmat_subrand(i,_)=arichBOTHmat_subrand(i+2*n_sampA,_);
}
//get unique alleles for both populations
for (int rich = 0; rich<nSNP;rich++){
NumericMatrix::Column inputBOTH_sub = arichBOTHmat_subrand(_,rich); 
NumericVector xBOTH_sub = clone<NumericVector>(inputBOTH_sub);
NumericVector tally_allelesBOTH_sub(2*n_sampA+2*n_sampB); //this is a vector that holds the unique alleles
int tallycounter_sub=1;

int nallsBOTH_sub=0;
typedef std::map<double,int> imap ;
typedef imap::value_type pair ;
imap index ;
int n = xBOTH_sub.size() ;
double current, previous = xBOTH_sub[0] ;
index.insert( pair( previous, 0 ) );
imap::iterator it = index.begin() ;
for( int i=1; i<n; i++){
current = xBOTH_sub[i] ;
if( current == previous ){
xBOTH_sub[i] = current + ( ++(it->second) / 100.0 ) ;
} else {
it = index.find(current) ;
if( it == index.end() ){
it = index.insert(
current > previous ? it : index.begin(),
pair( current, 0 )
) ;
} else {
xBOTH_sub[i] = current + ( ++(it->second) / 100.0 ) ;
}
previous = current ;
}

if (xBOTH_sub[i]-inputBOTH_sub[i]==0)
{nallsBOTH_sub=nallsBOTH_sub+1;
tally_allelesBOTH_sub[tallycounter_sub]=inputBOTH_sub[i];
tallycounter_sub=tallycounter_sub+1;
}
}
tally_allelesBOTH_sub[0]=inputBOTH_sub[0];

tallyBOTHmat_sub(_,rich)=tally_allelesBOTH_sub;
} //end of rich BOTH loop (over 10 microsats)

List tallist_sub(nSNP); //list of 10 that will hold allele frequencies for each of 10 microsatellites.
//First column is the actual allele size, second is freq from popA, third is freq popB
//now set up a table to get allele frequencies for popA and popB
for (int rich=0;rich<nSNP;rich++){
NumericMatrix::Column inputA_sub=arichAmat_subrand(_,rich);
NumericMatrix::Column inputB_sub=arichBmat_subrand(_,rich);
NumericMatrix::Column tal_sub=tallyBOTHmat_sub(_,rich); 

for(int i=0;i<tal_sub.size();i++){  //figure out how many unique alleles are at this microsat (from tal)
if(tal_sub[i]!=0){count_sub[rich]=count_sub[rich]+1;}
}

NumericVector Counter_sub=clone(count_sub);
int counter_sub=Counter_sub[rich];

NumericVector taltrunc_sub(counter_sub);  //will hold all unique alleles from this microsat (both pops)
NumericMatrix tallyAB_sub(counter_sub,3);  //matrix that has alleles, freq at popA, freq at popB
NumericVector howmanyA_sub(counter_sub);  //number of alleles for this population at each allele
NumericVector howmanyB_sub(counter_sub);


for(int i=0;i<counter_sub;i++){  //counter is the number of unique alleles at a locus
taltrunc_sub[i]=tal_sub[i];
int counterA_sub=0;  //a counter for number of unique alleles at each locus
int counterB_sub=0;

for (int j=0;j<2*n_sampA;j++){
if (inputA_sub[j]==taltrunc_sub[i])//go through all alleles to see how many match this unique one
{counterA_sub=counterA_sub+1;}
}
howmanyA_sub[i]=counterA_sub;


for (int j=0;j<2*n_sampB;j++){
if (inputB_sub[j]==taltrunc_sub[i])
{counterB_sub=counterB_sub+1;}
}
howmanyB_sub[i]=counterB_sub;

} //end of counter

tallyAB_sub(_,0)=taltrunc_sub;
tallyAB_sub(_,1)=howmanyA_sub/(2*n_sampA);
tallyAB_sub(_,2)=howmanyB_sub/(2*n_sampB);

arichsampA_ten[rich]=0;
arichsampB_ten[rich]=0;

for(int i=0;i<taltrunc_sub.size();i++){
if (howmanyA_sub[i]!=0){arichsampA_ten[rich]=arichsampA_ten[rich]+1;}
if (howmanyB_sub[i]!=0){arichsampB_ten[rich]=arichsampB_ten[rich]+1;}
}


tallist_sub[rich]=tallyAB_sub;
}//end of rich

arichsampA_vec[ngen+k]=std::accumulate(arichsampA_ten.begin(),arichsampA_ten.end(),0.0)/nSNP;
arichsampB_vec[ngen+k]=std::accumulate(arichsampB_ten.begin(),arichsampB_ten.end(),0.0)/nSNP;

//create talmat, which has unique alleles first column then freqs for popA and B
NumericMatrix talmat_sub(std::accumulate(count_sub.begin(),count_sub.end(),0.0),5);

for (int i=0;i<nSNP;i++){
int talcount_sub=0;
NumericMatrix taltmp_sub(count_sub[i],3);
taltmp_sub=as<SEXP>(tallist_sub[i]);
for (int j=std::accumulate(count_sub.begin(),count_sub.begin()+i,0.0);j<std::accumulate(count_sub.begin(),count_sub.begin()+i+1,0.0);j++)
{
talmat_sub(j,_)=taltmp_sub(talcount_sub,_);
talcount_sub=talcount_sub+1;
}
}
talmat_sub2=clone(talmat_sub);

double hetsA_sub;
double hetsB_sub;
for (int i=0;i<nSNP;i++){

for (int j=std::accumulate(count_sub.begin(),count_sub.begin()+i,0.0);j<std::accumulate(count_sub.begin(),count_sub.begin()+i+1,0.0);j++)
{
hetsA_sub=0;hetsB_sub=0;
for (int k=0;k<100;k++){
if ((arichAmat_sub((2*k),i)!=arichAmat_sub(((2*k)+1),i))&((arichAmat_sub((2*k),i)==talmat_sub2(j,0))|(arichAmat_sub(((2*k)+1),i)==talmat_sub2(j,0)))){hetsA_sub=hetsA_sub+1;}
if ((arichBmat_sub((2*k),i)!=arichBmat_sub(((2*k)+1),i))&((arichBmat_sub((2*k),i)==talmat_sub2(j,0))|(arichBmat_sub(((2*k)+1),i)==talmat_sub2(j,0)))){hetsB_sub=hetsB_sub+1;}
}
talmat_sub2(j,3)=hetsA_sub/100;
talmat_sub2(j,4)=hetsB_sub/100;
}}

//talmat: GET FST for SUBSAMPLE this function finds Fst between 2 populations, multiple loci and multiple alleles.
double n_bar_sub=0.5*(n_sampA+n_sampB); //sample size - can change this. already declared at 1746
double r_sub=2;
NumericVector p_bar_sub(talmat_sub2.nrow());
NumericMatrix::Column sampmatA_sub =talmat_sub2(_,1);
NumericMatrix::Column sampmatB_sub =talmat_sub2(_,2);
NumericVector s2_sub(talmat_sub2.nrow());
NumericVector h_bar_sub(talmat_sub2.nrow());

for (int i=0;i<talmat_sub2.nrow();i++){
h_bar_sub[i]=(talmat_sub2(i,3)+talmat_sub2(i,4))/2;
}

NumericVector ones_sub(talmat_sub2.nrow());
ones_sub=rep(1,talmat_sub2.nrow());
p_bar_sub=(n_sampA*sampmatA_sub+n_sampB*sampmatB_sub)/(2*n_bar_sub);  //each entry is average sample frequency of an allele
s2_sub=((n_sampB*(pow((sampmatB_sub-p_bar_sub),2)))+(n_sampA*pow((sampmatA_sub-p_bar_sub),2)))/n_bar_sub;//sample variance of allele freqs over pops (for each allele)

double nc_sub=((r_sub*n_bar_sub)-((pow(n_sampA,2)/(r_sub*n_bar_sub))+(pow(n_sampB,2)/(r_sub*n_bar_sub))))/(r_sub-1);  //same as n_bar

NumericVector a_sub=(n_bar_sub/nc_sub)*(s2_sub-((1/(n_bar_sub-1))*((p_bar_sub*(1-p_bar_sub))-((0.5)*s2_sub)-(0.25*h_bar_sub))));
NumericVector dvec_sub=((2*n_bar_sub)/((2*n_bar_sub)-1))*((p_bar_sub*(1-p_bar_sub))-(1/2)*s2_sub);
NumericVector b_sub=(n_bar_sub/(n_bar_sub-1))*((p_bar_sub*(1-p_bar_sub))-((1/2)*s2_sub)-(h_bar_sub*(2*n_bar_sub-1)/(4*n_bar_sub)));
NumericVector c_sub=h_bar_sub/2;
NumericVector aplusdvec_sub=a_sub+b_sub+c_sub;
double fst_sub=std::accumulate(a_sub.begin(),a_sub.end(),0.0)/std::accumulate(aplusdvec_sub.begin(),aplusdvec_sub.end(),0.0); //the better one
fst_subvec[sigfst]=fst_sub;

//if (sigfst==0&k==620|k==640|k==680){
//double C2_sub=r_sub*(1-(nc_sub/n_bar_sub));
//NumericVector ps_sub=p_bar_sub*(1-p_bar_sub);
//NumericVector theta_sub=(s2_sub-((1/((2*n_bar_sub)-1))*((p_bar_sub*(1-p_bar_sub))-s2_sub)))/(((1-((2*n_bar_sub*C2_sub)/(((2*n_bar_sub)-1)*r_sub)))*(p_bar_sub*(1-p_bar_sub)))+(1+((2*n_bar_sub*C2_sub)/((2*n_bar_sub-1)*r_sub)))*(s2_sub/r_sub));             

//NumericVector theta1_sub=(s2_sub-(1/(n_bar_sub-1)*((p_bar_sub*(1-p_bar_sub))-s2_sub-(h_bar_sub/4))))/((((1-((n_bar_sub*C2_sub)/((n_bar_sub-1)*r_sub)))*(p_bar_sub*(1-p_bar_sub)))+(((1+n_bar_sub*C2_sub)/((n_bar_sub-1)*r_sub))*(s2_sub/r_sub))+((C2_sub/(r_sub*(n_bar_sub-1)))*(h_bar_sub/4))));


//Wrights_fst_vec_sub1[ngen+k]=std::accumulate(theta1_sub.begin(),theta1_sub.end(),0.0)/talmat_sub2.nrow();
//Wrights_fst_vec_sub[ngen+k]=std::accumulate(theta_sub.begin(),theta_sub.end(),0.0)/talmat_sub2.nrow();
//Wrights_simple_sub[ngen+k]=std::accumulate(s2_sub.begin(),s2_sub.end(),0.0)/std::accumulate(ps_sub.begin(),ps_sub.end(),0.0);
//}

if(sigfst>0){
if(fst_sub>=fst_subvec[0]){fstsig[0]=fstsig[0]+1;}//count how many times random value is greater than true
}
}//end of sigfst

NumericVector outfstsig=clone(fstsig);

if (k==640&mgmt_counter==0)		{sig_submat60(zcounter,0)=0.01*outfstsig[0];//this gives significance of each Fst value
fst_submat60(zcounter,0)=fst_subvec[0];//this gives subsampled Fst
}

if (k==680&mgmt_counter==0)		{sig_submat99(zcounter,0)=0.01*outfstsig[0];//this gives significance of each Fst value
fst_submat99(zcounter,0)=fst_subvec[0];//this gives subsampled Fst
}

if (k==620&mgmt_counter==0)	{sig_submat9(zcounter,0)=0.01*outfstsig[0];
fst_submat9(zcounter,0)=fst_subvec[0];
}

//arichAmat_sub=clone(arichBmat_sub); //each time shift so that we do not double sample		
//n_sampA=n_sampB;

//}//end of nine comparisons
		
} //end of if nsamp>Ntotal something like that
}//end of k==9 or 60 and NtotalA and NtotalB>0

}//end of k= 1 to Nyrs popB 

out1(zcounter,_)=mean_arichA;
out2(zcounter,_)=mean_arichB;
out25(zcounter,_)=fst_vec;

} //end of mgmt_counter
} //end of zcounter
}//end if popdyn==1

//mgmt==0 (combined)
outmat[0]=out1; //mean_arichA;
outmat[1]=out2; //mean_arichB;
outmat[2]=out3; //R_l0A
outmat[3]=out4; //R_l0B
outmat[4]=out5; //NmigA_vec
outmat[5]=out6; //NmigB_vec
outmat[6]=out7; //TACA
outmat[7]=out8; //TACB
outmat[8]=out9; //FishAvec (fishing mortality calculated for population A - either sep or comb)
outmat[9]=out10; //FishBvec (fishing mortality calculated for population B - either sep or comb)

outmat[13]=out14; //TAC total catch for both pops combined mgmt only?
outmat[14]=out15;//NtotalA;
outmat[15]=out16;//NtotalB;
outmat[20]=out21; //holds SSB/initialSSB popA
outmat[21]=out22; //holds SSB/initial SSB popB
outmat[22]=out23; //holds estSSBA/SB40 popA
outmat[23]=out24; //holds estSSBB/SB40 popA
outmat[24]=out25;//fst_vec;
outmat[25]=Wrights_fst_vec;
outmat[26]=Wrights_fst_vec1;
outmat[27]=fst_vec2;
outmat[28]=Wrights_simple;
outmat[29]=ps;;
outmat[30]=Wrights_fst_vec_sub1;
outmat[31]=Wrights_simple_sub;
outmat[32]=s2;
//outmat[3]=allelesA;
//outmat[32]=allelesB;
outmat[33]=N_init_a1;
outmat[34]=N_init_a2;
outmat[35]=allA;
outmat[36]=allB;
outmat[37]=NeAvec;
outmat[38]=NeBvec;

//outmat[39]=popA;
outmat[40]=nmig_vec;

outmat[210]=fst_submat9; //this is the subsampled Fst, 1st column is true
outmat[211]=sig_submat9; //holds significances
outmat[212]=fst_submat60; //this is the subsampled Fst, 1st column is true
outmat[213]=sig_submat60; //holds significances
outmat[214]=fst_submat99; //this is the subsampled Fst, 1st column is true
outmat[215]=sig_submat99; //holds significances
outmat[216]=arichAmat_sub620;
outmat[217]=arichBmat_sub620;
outmat[218]=arichAmat_sub640;
outmat[219]=arichBmat_sub640;
outmat[220]=arichAmat_sub680;
outmat[221]=arichBmat_sub680;
outmat[230]=inputs;

return(outmat);

'
run2pops = cxxfunction(signature(INPUTS="numeric"), body = run8pops,plugin = "RcppArmadillo")
#out=run2pops(c(115500,0.2,110,100,2500,40,0.2389,1,1,1,3000,1,1000,1))

#SSBinit,pop1_prop,Nyrs,nsamp,ngen,nmig,fishmort,samptype,nsims,N,Yrs,mmt,dpow (power of distance),Rmat,ages, popdyn (1 for yes 0 for no),nSNP (# of markers)

#out=run2pops(c(5457,0.5,700,100,1,.34,0.12,1,1,3000,1000,1,1.15,1,7,1,13,85,100000,0))
#save(out,file="silverLL6_1.RData")

#out=run2pops(c(5457,0.5,700,100,1,.34,0.12,1,1,3000,1000,1,1.15,1,7,1,13,85,100000,0))
#save(out,file="silverLL6_2.RData")

#out=run2pops(c(5457,0.5,700,100,1,.34,0.12,1,1,3000,1000,1,1.15,1,7,1,13,85,100000,0))
#save(out,file="silverLL6_3.RData")

#out=run2pops(c(5457,0.5,700,100,1,.34,0.12,1,1,3000,1000,1,1.15,1,7,1,13,85,100000,0))
#save(out,file="silverLL6_4.RData")

#out=run2pops(c(5457,0.5,700,100,1,.34,0.12,1,1,3000,1000,1,1.15,1,7,1,13,85,100000,0))
#save(out,file="silverLL6_5.RData")


out=run2pops(c(ssba,0.5,700,100,1,.34,0.12,1,1,3000,1000,1,1.15,1,7,1,13,mig,100000,0))

save(out,file=paste("silverLLL_",ssba,"_",mig,"_",seed,".RData",sep="")) 

