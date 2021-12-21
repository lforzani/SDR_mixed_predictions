clear all;

p = 6; 
r = 2;
d = 1;
nreps = 100;
ns = [50 100 200 300 500];
nmax = max(ns);
N = length(ns);

mu = zeros(1,p);
I = eye(p);
alpha = zeros(p,1);
alpha((p/2):p) = 1/sqrt(p/2);
rho = 0.5;
Delta = (I+rho*alpha*alpha')*5;
        
n1 = 1;
n2 = 2;
epsi = [n1,n2];

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

 