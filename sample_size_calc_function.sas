/*******************************************************************************
    NAME    : sample_size_calc_function.sas
    TITLE   : SAS user-defined functions to calculate statistical power for cluster randomized trial with survival outcome
    PRODUCT : SAS R9.4
    AUTHOR  : Y.Inaba
    DATE    : 2017/09/07
*******************************************************************************/
***Proc FCMP;
libname fcn ".\Function";
libname data ".\data";
proc fcmp outlib = fcn.STD_Function.sample_size;

***Cluster-Level Weighted Log-Rank test;
**alpha:significant level(usually 0.05);
**n:number of clusters;
**D:the marginal probability of observing an event;
**p1:the limiting proportions of units assigned to treatments 1 and 2;
**gamma:local alternative(A marginal unweighted proportional hazards local alternative yields U = 1 for the unweighted log-rank
statistic);
**rho:intraclass correlation coeÆcient;
**m:common cluster size;
function sample_size_CLR(alpha,n,D,p1,gamma,rho,m);
Z_a=Quantile("NORMAL",1-(alpha/2));
p2=1-p1;
a1=p1;
a2=p2;
_q=sqrt((n*m*D*a1*a2*(gamma**2))/(1+(m-1)*rho))-Z_a;
beta=CDF("NORMAL",_q);
return(beta);
endsub;


**SharedFrailtyModel;
function sample_size_SFM(alpha,N,K,PC,beta1,CV);
Z_a=Quantile("NORMAL",1-(alpha/2));*print Z_a;
_q=sqrt((N * K * (1-PC) * (beta1**2) * ((1-exp(beta1))**2) )
                     / ((2 * ((1-exp(beta1))**2)) + (K * (beta1**2) * (1-PC) * (CV**2) * (1+exp(2*beta1)) ))
               )-Z_a;
beta=CDF("NORMAL",_q);
return(beta);
endsub;

**LWA model;
function sample_size_LWA(n,alpha,K,phi,P1,P2,theta,rho);
    Z_a=Quantile("NORMAL",1-(alpha/2),0,1);
    d=((K*((1-P1)+phi*(1-P2))*n)/(1+phi));
    _q=sqrt( ((d*phi*((1-theta)**2)))/(((1+phi*theta)**2)*(1+(K-1)*rho)) )-Z_a;
    beta=CDF("NORMAL",_q);
    return(beta);
endsub;

run;


***Example;

options cmplib=fcn.STD_Function;

data test;
infile datalines;
length n    alpha   K   phi P1  P2  theta   rho PC  beta1   CV  D   p1  gamma 8;
input n alpha   K   phi P1  P2  theta   rho PC  beta1   CV  D   p1  gamma ;
datalines;
10  0.05    10  1   0.2 0.24    0.64    0.1 0.22    0.60653066  0.11    0.78    0.5 0.5
10  0.05    30  1   0.2 0.24    0.64    0.1 0.22    0.60653066  0.11    0.78    0.5 0.5
10  0.05    50  1   0.2 0.24    0.64    0.1 0.22    0.60653066  0.11    0.78    0.5 0.5
30  0.05    10  1   0.2 0.24    0.64    0.1 0.22    0.60653066  0.11    0.78    0.5 0.5
30  0.05    30  1   0.2 0.24    0.64    0.1 0.22    0.60653066  0.11    0.78    0.5 0.5
30  0.05    50  1   0.2 0.24    0.64    0.1 0.22    0.60653066  0.11    0.78    0.5 0.5
50  0.05    10  1   0.2 0.24    0.64    0.1 0.22    0.60653066  0.11    0.78    0.5 0.5
50  0.05    30  1   0.2 0.24    0.64    0.1 0.22    0.60653066  0.11    0.78    0.5 0.5
50  0.05    50  1   0.2 0.24    0.64    0.1 0.22    0.60653066  0.11    0.78    0.5 0.5
10  0.05    10  1   0.2 0.24    0.64    0.5 0.22    0.60653066  0.11    0.78    0.5 0.5
10  0.05    30  1   0.2 0.24    0.64    0.5 0.22    0.60653066  0.11    0.78    0.5 0.5
10  0.05    50  1   0.2 0.24    0.64    0.5 0.22    0.60653066  0.11    0.78    0.5 0.5
30  0.05    10  1   0.2 0.24    0.64    0.5 0.22    0.60653066  0.11    0.78    0.5 0.5
30  0.05    30  1   0.2 0.24    0.64    0.5 0.22    0.60653066  0.11    0.78    0.5 0.5
30  0.05    50  1   0.2 0.24    0.64    0.5 0.22    0.60653066  0.11    0.78    0.5 0.5
50  0.05    10  1   0.2 0.24    0.64    0.5 0.22    0.60653066  0.11    0.78    0.5 0.5
50  0.05    30  1   0.2 0.24    0.64    0.5 0.22    0.60653066  0.11    0.78    0.5 0.5
50  0.05    50  1   0.2 0.24    0.64    0.5 0.22    0.60653066  0.11    0.78    0.5 0.5
10  0.05    10  1   0.2 0.24    0.64    0.9 0.22    0.60653066  0.11    0.78    0.5 0.5
10  0.05    30  1   0.2 0.24    0.64    0.9 0.22    0.60653066  0.11    0.78    0.5 0.5
10  0.05    50  1   0.2 0.24    0.64    0.9 0.22    0.60653066  0.11    0.78    0.5 0.5
30  0.05    10  1   0.2 0.24    0.64    0.9 0.22    0.60653066  0.11    0.78    0.5 0.5
30  0.05    30  1   0.2 0.24    0.64    0.9 0.22    0.60653066  0.11    0.78    0.5 0.5
30  0.05    50  1   0.2 0.24    0.64    0.9 0.22    0.60653066  0.11    0.78    0.5 0.5
50  0.05    10  1   0.2 0.24    0.64    0.9 0.22    0.60653066  0.11    0.78    0.5 0.5
50  0.05    30  1   0.2 0.24    0.64    0.9 0.22    0.60653066  0.11    0.78    0.5 0.5
50  0.05    50  1   0.2 0.24    0.64    0.9 0.22    0.60653066  0.11    0.78    0.5 0.5
;
run;
data power;*(keep=N K rho comp:);
    set test;
    comp_power_CLR=round(sample_size_CLR(alpha,n,D,p1,log(theta),rho,K)*100,0.1);
    comp_power_LWA=round(sample_size_LWA(n*2,alpha,K,phi,P1,P2,theta,rho)*100,0.1);
    comp_power_SFM=round(sample_size_SFM(alpha,N/2,K,PC,beta1,CV)*100,0.1);
run;
