clear all;
%continua binaria estimacion d=1
n =2000;
q = 10; %numero de discretas
p=20; %numeros de continuas
r = 5;
d=1;
nreps = 100;

k=[100 200 300 500 750];

N = length(k);
nmax = max(k);

d1=1;
d2=1;


mu = zeros(1,p);
I = eye(p);
alpha=zeros(p,1);
alpha(floor(p/2):p) = 1/sqrt(p/2);
alpha=alpha';
rho = 0.55;
alpha=alpha';
Delta = (I+rho*alpha*alpha')*5;
epsi=ones(r,1);
epsi=epsi';
A=Delta*alpha*epsi;
 


Gamma1 = zeros(q,q);
Gamma1(1,1) = 1;
Gamma1(2,2) = 1;
Gamma1(3,3) = 1;
Gamma1(4,4) = 1;
Gamma1(1,2) = 30;
Gamma1(2,1) = 30;
Gamma1(1,3) = 5;
Gamma1(3,1) = 5;
Gamma1(3,4) = 30;
Gamma1(4,3) = 30;
Gamma1(2,3) = 10;
Gamma1(3,2) = 10;
for i=2:q
    Gamma1(i,i)=1;
    Gamma1(i-1,i)=30;
    Gamma1(i,i-1)=30;
end

Gamma1=Gamma1/sqrt(sum(sum(Gamma1)));
Gamma=2*[Gamma1,Gamma1,Gamma1,Gamma1,Gamma1];

Gamma=3*[Gamma1,Gamma1,Gamma1,Gamma1,Gamma1];

Gannas=reshape(Gamma,[q,q,r]);



 
beta=zeros(p,q); 
mubeta=zeros(1,q);
Ibeta=eye(q); 
for ii = 1:p
  beta(ii,:) = 1/10;
end
beta(:,1:(q/2))=0;
beta=beta;

 
 

auxx2 = ones(q,q);
True_suboptimal = [alpha ; -beta'*alpha; diag(Gamma1); Gamma1(find(tril(auxx2,-1)))];

True_optimal=[alpha;diag(Gamma1)-beta'*alpha;Gamma1(find(tril(auxx2,-1)))];

 
 

%%
g = r+1;
probs = ones(g,1)/g;
YY1 =  randsample(g, n,true, probs);
fycent1=get_fyZ(YY1);
[HH1 , XX1] = GenData_NonNormal(Gannas,beta,Delta,A,fycent1,'mixtos');
[T_optimal,T_suboptimal] = getT_optimal(XX1,HH1);
%%

       Proy_True_optimal = T_optimal*True_optimal;
       True_suboptimal_cont=True_suboptimal(1:(p+q));
       True_suboptimal_disc=True_suboptimal((p+q+1):((p+q+q*(q+1)/2)));
       T_suboptimal_cont=T_suboptimal(:,(1:(p+q)));
       T_suboptimal_disc=T_suboptimal(:,(p+q+1):((p+q+q*(q+1)/2)));
       Proy_True_suboptimal_cont = T_suboptimal_cont*True_suboptimal_cont;
       Proy_True_suboptimal_disc = T_suboptimal_disc*True_suboptimal_disc

for jj = 1:nreps
    jj
    YY =  randsample(g, n,true, probs);
    fycent=get_fyZ(YY);
    [HH , XX] =  GenData_NonNormal(Gannas,beta,Delta,A,fycent,'mixtos');
    [TTT_optimal,~] = getT_optimal(XX,HH);

    for mm = 1:length(k)
       Y=YY(1:k(mm));
       H=HH(1:k(mm),:);
       X=XX(1:k(mm),:);    
       TT_optimal =TTT_optimal(1:k(mm),:);
      
      %[redu_optimal,redu_suboptimal,fycent,~,~,Tauhatraro0, Tauhatraro, Ahat, betahat, Deltahat] = EM4mixture_continua_binaria(X,H,Y,'disc');
      
       %We are proving unsupervised with our method
       YYY = ones(size(Y)); 
       %[redu_optimal_us,redu_suboptimal_us,fycent_us,~,~,Tauhatraro0_us, Tauhatraro_us, Ahat_us, betahat_us, Deltahat_us] = EM4mixture_continua_binaria(X,H,YYY,1);
       
       
      
        datos=[X,H];
 
       
        %%%%%%aotuomatic choosen of a with true d
         % suboptimal
         [red_est_suboptimal1,red_est_suboptimal2,betahat] = rr4mix_suboptimal(X,H,Y,d1,d2,'disc','auto');
          
         
         %optimal choose a with true d
         [red_est_optimal] = rr4mix_optimal(X,H,Y,d,'disc','auto');
  
         
        
         Proy_est_optimal = T_optimal*red_est_optimal;
          
         Proy_est_suboptimal11_cont = T_suboptimal_cont(:,1:p)*red_est_suboptimal1(1:p)-T_suboptimal_cont(:,(p+1):(p+q))*betahat'*red_est_suboptimal1(1:p);
         Proy_est_suboptimal11_disc = T_suboptimal_disc*red_est_suboptimal2;
      

         %%Computing PCA and PFC on T_optimal (i.e. over X,H and H_iH_j)
         [coefPCA, scoresPCA] = pca(TT_optimal);
         reduPCA=coefPCA(:,1:d);
         Proy_PCA = T_optimal*reduPCA; 
         
         %PFC taking all variables as continuous  
       [Proy_est_PFCTT,redu_est_PFC] = ldr(Y,TT_optimal,'pfc','disc',d);
        Proy_est_PFC = T_optimal*redu_est_PFC;

           %%angulos suboptimals
         
              [RDsubop1auto_F(jj,mm), RDsubop1auto_2(jj,mm)]=diferencia(red_est_suboptimal1,True_suboptimal_cont);
       
              
              [RDsubop2auto_F(jj,mm), RDsubop2auto_2(jj,mm)] = diferencia(red_est_suboptimal2,True_suboptimal_disc);
      
              
              
%%%%optimal
         
              
             [RDoptimalauto_F(jj,mm), RDoptimalauto_2(jj,mm)] = diferencia(red_est_optimal,True_optimal);
              
         
             [RD_PCA_F(jj,mm), RD_PCA_2(jj,mm)] = diferencia(reduPCA,True_optimal);
             
         
              
             [RD_PFC_F(jj,mm), RD_PFC_2(jj,mm)] = diferencia(redu_est_PFC,True_optimal);
             
              
             
              % Proyecciones
                % optimal
                   
                   [PDoptimalauto_F(jj,mm), PDoptimalauto_2(jj,mm)]  = diferencia(Proy_est_optimal,Proy_True_optimal);
                   
        
                   
         [PD_PCA_F(jj,mm), PD_PCA_2(jj,mm)]  = diferencia(Proy_PCA,Proy_True_optimal);

                   
         [PD_PFC_F(jj,mm), PD_PFC_2(jj,mm)]  = diferencia(Proy_est_PFC,Proy_True_optimal);

     
                 % suboptimal
              
              [PDsuboptimalauto_cont_F(jj,mm), PDsuboptimalauto_cont_2(jj,mm)]  = diferencia(Proy_est_suboptimal11_cont,Proy_True_suboptimal_cont);
              
              
              
              [PDsuboptimalauto_disc_F(jj,mm), PDsuboptimalauto_disc_2(jj,mm)] = diferencia(Proy_est_suboptimal11_disc,Proy_True_suboptimal_disc);
       
     
        
     
     end

end

              
              
