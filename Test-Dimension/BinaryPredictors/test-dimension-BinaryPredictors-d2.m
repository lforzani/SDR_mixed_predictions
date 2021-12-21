clear all;

q = 5; 
r = 2;
d = 2;
nreps = 500;
ns = 1000;
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

for i=2:q
    Gamma1(i,i) = 1;
    Gamma1(i-1,i) = 30;
    Gamma1(i,i-1) = 30;
end

Gamma1 = Gamma1/sqrt(sum(sum(Gamma1)));

 
Gamma2=zeros(q,q);
 for i=1:q
    Gamma2(i,i) = 1;
end


Gamma = [Gamma1, Gamma2];
Gannas = reshape(Gamma,[q,q,r]);

g = r+1;

probs = ones(g,1)/g;

dtest1 = zeros(nreps,1); 
dtest2 = zeros(nreps,1);

D1 = zeros(nreps, N);
D2 = zeros(nreps, N);
    
for j = 1:nreps

    YY =  randsample(g, nmax, true, probs);

    fycent = get_fyZ(YY);
    
              [HH] = GenDataBinariasOPv2(nmax, fycent, Gannas); 

    for mm = 1:N
        
        n = ns(mm);

        Y = YY(1:n);
        
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
