clear all;

rng(1000)
 
% set one
setone = load('setone.csv');
p = 7;
q = 3;
datos = [setone(:,[1 2 3 4 5 7]) log(setone(:,6)) setone(:,8:11)];

% set two
%settwo = load('settwo.csv');
%settwo = [settwo(:,1:3), log(settwo(:,4)), settwo(:,5:8)];
%p = 4;
%q = 3;
%datos = settwo;


% set three
setthree = load('setthree.csv');
p = 8;
q = 4;
datos = [setthree(:,[2 7 8]), log(setthree(:,[1 3 4 5 6])), setthree(:,9:13)];

% set four
setfour = load('setfour.csv');
p = 6;
q = 5;
datos = [setfour(:,1:4), log(setfour(:,5:6)), setfour(:,7:12)];

n = size(datos,1);

X = datos(:,1:p);
H = datos(:,(p+1):(p+q));
Y = datos(:,(p+q+1));

% for the setfour
% H = datos(:,(p+1):(p+3));

type = 'disc';
model = 'logit';

[MSE_optimal, sd_optimal, MSE_suboptimal, sd_suboptimal, MSE_all, sd_all, MSE_PFC, sd_PFC, ProyPFC_tr,ProyPFC_ts, ProyOptimal_tr, ProyOptimal_ts,ProySuboptimal_tr1, ProySuboptimal_ts1,ProySuboptimal_tr2, ProySuboptimal_ts2, lambdas_optimal]= cvPFCmix_LOOCV(Y,X,H,type,model);
lambdaopt=mean(lambdas_optimal);


%Compute AUC using overall sample
n = size(Y,1);
p = size(X,2);
q = size(H,2);
kk = q*(q+1)/2;
dim = 1;
dim1 = 1;
dim2 = 1; 
red_est_optimal = rr4mix_optimalLR(X,H,Y,dim,type,lambdaopt);
[red_est_suboptimal1,red_est_suboptimal2] = rr4mix_suboptimal(X,H,Y,dim1,dim2,type,0);   

%Compute projections
Hraro = zeros(n,kk);
for j = 1:n
auxts = H(j,:)'*H(j,:);
Hraro(j,:) = [diag(auxts);auxts(find(tril(ones(q,q),-1)))];
end
T_optimal = [X,Hraro];
T_suboptimal1 = [X,H];
T_suboptimal2 = Hraro;

V = [X H];

Proy_est_optimal= T_optimal*red_est_optimal;
Proy_est_suboptimal1 = T_suboptimal1*red_est_suboptimal1; 

if(size(red_est_suboptimal2,2) == 0)
Proy_est_suboptimal = Proy_est_suboptimal1;
else     
Proy_est_suboptimal2 = T_suboptimal2*red_est_suboptimal2; 

Proy_est_suboptimal = [Proy_est_suboptimal1,Proy_est_suboptimal2];
end

%fit logit models
m_optimal = glmfit(Proy_est_optimal,Y,'binomial');
m_suboptimal = glmfit(Proy_est_suboptimal,Y,'binomial');
mm = glmfit(V, Y,'binomial');

hatprobas_optimal= glmval(m_optimal, Proy_est_optimal, 'logit');
[~,~,~,AUCoptimal] = perfcurve(Y,hatprobas_optimal,1);

hatprobas_suboptimal=glmval(m_suboptimal, Proy_est_suboptimal, 'logit');
[~,~,~,AUCsuboptimal] = perfcurve(Y,hatprobas_suboptimal,1);

hatprobas_full=glmval(mm, V, 'logit');
[~,~,~,AUCfull] = perfcurve(Y,hatprobas_full,1);

 
[MSE_optimal,MSE_suboptimal,MSE_all]

red_est_optimal
%red_est_suboptimal1 
%red_est_suboptimal2
 
[AUCoptimal,AUCsuboptimal,AUCfull]
    