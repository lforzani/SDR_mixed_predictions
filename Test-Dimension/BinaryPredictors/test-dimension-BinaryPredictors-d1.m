clear all;

q = 8; 
r = 2;
d = 1;
nreps = 500;
ns = [50 100 200 300 500];
N = length(ns);
nmax = max(ns);
d1 = 1;
d2 = 1;

n1 = 1;
n2 = 2;

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
 
 for i = 2:q
     Gamma1(i,i)=1;
     Gamma1(i-1,i)=30;
     Gamma1(i,i-1)=30;
 end
 
Gamma1 = Gamma1/sqrt(2004);
Gamma2 = 2*Gamma1;

Gamma = [Gamma1, Gamma2];
Gannas = reshape(Gamma,[q,q,r]);

auxx2 = ones(q,q);
True_suboptimal = [diag(Gamma1); Gamma1(find(tril(auxx2,-1)))];
True_optimal = [diag(Gamma1); Gamma1(find(tril(auxx2,-1)))];

g = r+1;

probs = ones(g,1)/g;

dtest1 = zeros(nreps,1); 
dtest2 = zeros(nreps,1);

D1 = zeros(nreps, N);
D2 = zeros(nreps, N);
    
for j = 1:nreps

    YY =  randsample(g, nmax, true, probs);

    fycent = get_fyZ(YY);
      
    for mm = 1:N
        
        n = ns(mm);

        Y = YY(1:n);
        
        [HH] = GenDataBinariasOPv2(nmax, fycent, Gannas); 

        H = HH(1:n,:);
    
        [T, tau0,redu,proj,fycent] = EM4mixture_binaria(H,Y,'disc');
        Vrcl = AsymptoticVariancebinarias(tau0,redu,fycent,q);

        [aux1,aux2] = testChid(n,Vrcl,redu, 0.05);

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
