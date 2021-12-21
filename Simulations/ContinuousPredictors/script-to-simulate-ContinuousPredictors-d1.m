clear all;
%Simulation for continuous predictor with d=1 
n = 2000;
p=20; 
r = 5;
d=1;
nreps = 100;
k=[100 200 300 500 750];
 
N = length(k);
nmax = max(k);

mu = zeros(1,p);
I = eye(p);
alpha=zeros(p,1);
alpha((1+floor(p/2)):p) = 1/sqrt(p/2);
alpha=alpha';
rho = 0.55;
alpha=alpha';

Delta=rho*ones(p,p);
for j=1:p
    Delta(j,j)=1;
end

Delta=Delta*5
 
Delta = (I+rho*alpha*alpha')*5;
AA = randn( p, p );  % random iid ~N(0,1)
oAA = orth( AA.' ).'; % orthogonal rows
nA = bsxfun( @rdivide, oAA, sqrt( sum( oAA.^2, 2 ) ) ); % normalize to unit length

Delta=nA*Delta*nA';
epsi=ones(r,1);
epsi=epsi';
A=Delta*alpha*epsi;
 
True_optimal=[alpha];

%To count true selected variables
%%True zeros
 TZ = find(alpha==0);
%%True selected variables (non zeros)
 TNZ = find(alpha);  
%Vectors to count coincidences and comput TPR and TNR (Table 2) 
CoincidZ=zeros(nreps,length(k));
CoincidNZ =zeros(nreps, length(k));
 
%
g = r+1;
probs = ones(g,1)/g;
YY1 =  randsample(g, n,true, probs);
fycent1=get_fyZ(YY1);
[XX1] = GenDataContinuouss(fycent1,Delta,A);
         Proy_True_optimal = XX1*True_optimal;
         
for jj = 1:nreps
    jj
    YY =  randsample(g, n,true, probs);
    fycent=get_fyZ(YY);
    [XX] = GenDataContinuouss(fycent,Delta,A);
    for mm = 1:length(k)
       Y=YY(1:k(mm));
       
       X=XX(1:k(mm),:);
       
 %  Parameter estimation without reduction       
 [T_optimal,redu_optimal,proj,Ahat,Deltahat,fycent] = EM4mixture_Continuous(X,Y,'disc');

         
 % Optimal reduction with automatic variable selection
        
 [red_est_optimal] = rr4cont(X,Y,d,'disc','auto');
 
 % Compute the projecton        
 Proy_est_optimal = XX1*red_est_optimal;
           
 
 TrueZeros=sum(red_est_optimal(TZ)==0)/length(red_est_optimal(TZ));
 FalseZeros=1-(sum(red_est_optimal(TNZ)==0))/length(red_est_optimal(TNZ));
 CoincidZ(jj,mm) =TrueZeros;
 CoincidNZ(jj,mm) =FalseZeros;
         
 % Computing PCA over X to compare with unsupervised reduction method
         
 [coefPCA, scoresPCA] = pca(X);
 reduPCA=coefPCA(:,d);
 Proy_PCA = XX1*reduPCA; 

                 
 % Compute distante (Frobenious and Eucliean) between true and estimated reduction 
  
[Roptimalauto_F(jj,mm) Roptimalauto_2(jj,mm)] = diferencia(red_est_optimal,True_optimal);
        
[RD_PCA_F(jj,mm) RD_PCA_2(jj,mm)] = diferencia(reduPCA,True_optimal);
 
 % Compute distante (Frobenious and Eucliean) between out-of-sample projections
          
 [PDoptimalauto_F(jj,mm) PDoptimalauto_2(jj,mm)]= diferencia(Proy_est_optimal,Proy_True_optimal);
                             
 [PD_PCA_F(jj,mm) PD_PCA_2(jj,mm)]= diferencia(Proy_PCA,Proy_True_optimal);
        
    end

end


cd '~/Figures/ContinuousPredictors'


csvwrite('Continuous-d1-estimator.csv',Roptimalauto_2);
csvwrite('Continuous-d1-proyection.csv',PDoptimalauto_2);
          
csvwrite('Continuous-d1-estimator_PCA.csv',RD_PCA_2);
csvwrite('Continuous-d1-proyection_PCA.csv',PD_PCA_2);
          

