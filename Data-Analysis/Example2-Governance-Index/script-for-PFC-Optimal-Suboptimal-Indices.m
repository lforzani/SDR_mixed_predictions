%This script is to run the example of Governance index annotation.
% First we compute CGindexes using PCA, PFC and Optimal SDR, save them to
% work in R.
%Second, we run leave-on-out cross validation (using the function cvPFCmix_LOOOCV, in the internal folder) for PFC, Optimal  SDR and
%Subotimal SDR, exporting te results to run in R together with PCA and PCAmix.
%Nonparametric fits have be done using 'np' package of R.


clear all;

%set seed
rng(1234);


%Open-Read data set
datos = readtable('~/Data-Analysis/Example2-Governance-Index/GovernanceDataExtended2.csv', 'Delimiter', ';', 'HeaderLines',1);


% Specify the number of continuous (p) , binary predictors (q) y r 
p = 6;
q=11;
r=1;
X = table2array(datos(:,2:7));
H = table2array(datos(:,10:(10+(q-1))));
Y = table2array(datos(:,1));


%We run leave-one-out for PFC, Optimal and Suboptimal and save results of
%Projections to use in R

type = 1;
model = 'linear';

%[MSE_PFC, sd_PFC, MSE_optimal, sd_optimal, MSE_suboptimal, sd_suboptimal, MSE_all, sd_all] = cvPFCmix_LOOCV(Y,X,H,type,model);
[ MSE_optimal, sd_optimal, MSE_suboptimal, sd_suboptimal, MSE_all, sd_all, MSE_PFC, sd_PFC,ProyPFC_tr,ProyPFC_ts, ProyOptimal_tr, ProyOptimal_ts,ProySuboptimal_tr1, ProySuboptimal_ts1,ProySuboptimal_tr2, ProySuboptimal_ts2] = cvPFCmix_LOOCV(Y,X,H,type,model);





%Results of LOOCV of PFC, SDROptimal and SDRSubOptimal 

[MSE_PFC,MSE_optimal,MSE_suboptimal]
[sd_PFC,sd_optimal,sd_suboptimal]

%save projections
ProyPFC=[ProyPFC_ts; ProyPFC_tr];
ProyOptimal=[ProyOptimal_ts; ProyOptimal_tr];
ProySuboptimal1=[ProySuboptimal_ts1; ProySuboptimal_tr1];
ProySuboptimal2=[ProySuboptimal_ts2; ProySuboptimal_tr2];


csvwrite('ProyPFC.csv',ProyPFC);
csvwrite('ProyOptimal.csv',ProyOptimal);
csvwrite('ProySuboptimal1.csv',ProySuboptimal1);
csvwrite('ProySuboptimal2.csv',ProySuboptimal2);



