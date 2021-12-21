clear all;

p = 6;
r = 2;
d = 2;
nreps = 500;
ns = [50 100 200 300 500];
nmax = max(ns);
N = length(ns);

mu = zeros(1,p);
I = eye(p);

alpha1 = ones(p,1)/sqrt(p);
rho1 = 0.55;
alpha2 = ones(p,1)/sqrt(p);
rho2 = 0.25;
alpha1(p/2:p) = 0;
alpha2(p/2:p) = 0;

alpha2(1:floor(p/4)) =-1/sqrt(p);

alpha = [alpha1  alpha2]; 
alpha = orth(alpha);

Delta = (I+rho1*alpha(:,1)*alpha(:,1)'+rho2*alpha(:,2)*alpha(:,2)')*5;

n1 = mvnrnd(0,1);
n2 = mvnrnd(0,1);

n11 = mvnrnd(0,1);
n22 = mvnrnd(0,1);

epsi = [n1, n11; n2, n22]';
A = Delta*alpha*epsi;

g = r+1;

probs = ones(g,1)/g;

dtest1 = zeros(nreps,1); 
dtest2 = zeros(nreps,1);

D1 = zeros(nreps, N);
D2 = zeros(nreps, N);
    
for j = 1:nreps
    
    YY =  randsample(g, nmax, true, probs);

    fycent = get_fyZ(YY);
 
    XX = GenDataContinuas(fycent,Delta,A);

    for mm = 1:N
        
        n = ns(mm);

        Y = YY(1:n);
        X = XX(1:n,:);
        
        [T,redu_optimal,proj,Ahat,Deltahat,fycent] = EM4mixture_continua(X,Y,'disc');
        
        Vrcl = AsymptoticVariancecontinuas(Deltahat,Ahat,fycent);

        [aux1,aux2] = testChid(n,Vrcl,redu_optimal, 0.05);

        dtest1(j) = aux1;
        dtest2(j) = aux2;
   
    if (dtest1(j) == d)
        D1(j,mm) = 1;
    end
    
    if (dtest2(j) == d)
        D2(j,mm) = 1;
    end
    
    end
end

mean(D1, 1)
mean(D2, 1)
