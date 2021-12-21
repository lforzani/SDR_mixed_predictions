clear all;
%Simulation for non-normal continuous predictor with d=2 for Table 6 in
%Appendix D
n = 2000;
rng(1000);
p=20; 
r = 5;
d=2;
nreps = 100;
k=[100 200 300 500 750];

N = length(k);
nmax = max(k);

mu = zeros(1,p);
I = eye(p);
alpha1=zeros(p,1);
alpha1((1+floor(p/2)):p) = 1;
alpha2=alpha1;
alpha2(floor(3/4*p):p)=-1;
alpha=[alpha1, alpha2];
alpha=orth(alpha);
rho1 = 0.55;
rho2=.25;
 
Delta = (I+rho1*alpha(:,1)*alpha(:,1)'+rho2*alpha(:,2)*alpha(:,2)')*5;

AA = randn( p, p );  % random iid ~N(0,1)
oAA = orth( AA.' ).'; % orthogonal rows
nA = bsxfun( @rdivide, oAA, sqrt( sum( oAA.^2, 2 ) ) ); % normalize to unit length

Delta=nA*Delta*nA';

epsi=ones(r,d);
epsi(1:3,2)=0;
epsi=epsi';
A=Delta*alpha*epsi;
 
True_optimal=[alpha];

 %To count true selected variables
%%True zeros
TZ = find(alpha1==0 & alpha2==0);
%%selected variables (non zeros)
 TNZ = find(alpha1~=0 | alpha2~=0); 

%Vectors to count coincidences and comput TPR and TNR (Table 2) 
CoincidZ=zeros(nreps,length(k));
CoincidNZ =zeros(nreps, length(k));
 
%
g = r+1;
probs = ones(g,1)/g;
YY1 =  randsample(g, n,true, probs);
fycent1=get_fyZ(YY1);
[XX1] = GenDataContinuas_NonNormal(fycent1,Delta,A);
Proy_True_optimal = XX1*True_optimal;

for jj = 1:nreps
    jj
    YY =  randsample(g, n,true, probs);
    fycent=get_fyZ(YY);
    [XX] = GenDataContinuas_NonNormal(fycent,Delta,A);
    for mm = 1:length(k)
       Y=YY(1:k(mm));
       
       X=XX(1:k(mm),:);
       
 %  Parameter estimation without reduction       
 [T_optimal,redu_optimal,proj,Ahat,Deltahat,fycent] = EM4mixture_continua(X,Y,'disc');
     
 % Optimal reduction with automatic variable selection

  [red_est_optimal,~,~,lambda] = rr4cont(X,Y,d,'disc','auto');
  
   % Compute the projecton        
  Proy_est_optimal = XX1*red_est_optimal;
            
  TrueZeros=sum(red_est_optimal(TZ,1)==0 &red_est_optimal(TZ,2)==0)/length(red_est_optimal(TZ,2));
  FalseZeros=1-(sum(red_est_optimal(TNZ,1)==0 & red_est_optimal(TNZ,2)==0))/length(red_est_optimal(TNZ));
  CoincidZ(jj,mm) =TrueZeros;
  CoincidNZ(jj,mm) =FalseZeros;
  
  % Computing PCA over X to compare with unsupervised reduction method
         
  [coefPCA, scoresPCA] = pca(X);
  reduPCA=coefPCA(:,1:d);
  Proy_PCA = XX1*reduPCA; 
    
         
 % Compute distante (Frobenious and Eucliean) between true and estimated reduction 

 [Roptimalauto_F(jj,mm) Roptimalauto_2(jj,mm)]= diferencia(red_est_optimal,True_optimal);

 [RD_PCA_F(jj,mm) RD_PCA_2(jj,mm)] = diferencia(reduPCA,True_optimal);

 % Compute distante (Frobenious and Eucliean) between out-of-sample projections

 [PDoptimalauto_F(jj,mm) PDoptimalauto_2(jj,mm)]= diferencia(Proy_est_optimal,Proy_True_optimal);

 [PD_PCA_F(jj,mm) PD_PCA_2(jj,mm)]= diferencia(Proy_PCA,Proy_True_optimal);

        
    end

end





