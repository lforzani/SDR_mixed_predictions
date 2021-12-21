clear all;

q = 5; 
p = 6; 
r = 2;

nreps = 10;
ns = 100;
N = length(ns);
nmax = max(ns);
d1 = 1;
d2 = 1;

mu = zeros(1,p);
I = eye(p);
alpha = zeros(p,1);
alpha((p/2):p) = 1/sqrt(p/2);
rho = 0.55;
Delta = (I+rho*alpha*alpha')*5;
         
n1 = 1;
n2 = 2;
epsi = [n1,n2];

A = Delta*alpha*epsi;

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
Gamma2 = n2/n1*Gamma1;
Gamma = [Gamma1,Gamma2];

Gannas = reshape(Gamma,[q,q,r]);

%Generate beta
beta = zeros(p,q); 

mubeta = zeros(1,q);

Ibeta = eye(q); 
for ii = 1:p
    beta(ii,:) = 1/30;
end

auxx2 = ones(q,q);
True_suboptimal1 =[alpha ; -beta'*alpha];
True_suboptimal2 = [diag(Gamma1); Gamma1(find(tril(auxx2,-1)))];
True_optimal = [alpha; diag(Gamma1)-beta'*alpha; Gamma1(find(tril(auxx2,-1)))];

g = r+1;

probs = ones(g,1)/g;

dtestfirst1 = zeros(nreps,1);
dtestsecond1 = zeros(nreps,1);
dtestfirst2 = zeros(nreps,1);
dtestsecond2 = zeros(nreps,1);

Dfirst1 = zeros(nreps, N);
Dfirst2 = zeros(nreps, N);
Dsecond1 = zeros(nreps, N);
Dsecond2 = zeros(nreps, N);

for j = 1:nreps

    YY =  randsample(g, nmax, true, probs);

    fycent = get_fyZ(YY);

    [HH, XX] = GenDataBinariasContinuas(Gannas,fycent,Delta,A,beta);

    for mm = 1:N
        
        n = ns(mm);

        Y = YY(1:n);
        H = HH(1:n,:);
        X = XX(1:n,:);

        [redu_optimal,redu_suboptimal, fycent,T_optimal,T_suboptimal, Tauhatraro0, Tauhatraro, Ahat, betahat, Deltahat] = EM4mixture_continua_binaria(X,H,Y,'disc');

        Vrclfirst = AsymptoticVariancecontinuas(Deltahat,Ahat,fycent);

        redu_suboptimalfirst =redu_suboptimal(1:p,:);
        [auxfirst1,auxfirst2] = testChid(n,Vrclfirst,redu_suboptimalfirst, 0.05);

        s= size(redu_suboptimal);
        redu_suboptimalsecond =redu_suboptimal(p+q+1:s,:);

        Vrclsecond = AsymptoticVariancebinarias(Tauhatraro0,Tauhatraro,fycent,q);

        [auxsecond1,auxsecond2] = testChid(n,Vrclsecond,redu_suboptimalsecond, 0.05);

        dtestfirst1(j) = auxfirst1;
        dtestfirst2(j) = auxfirst2;

        dtestsecond1(j) = auxsecond1;
        dtestsecond2(j) = auxsecond2;
        
    if (dtestfirst1(j) == d1)
        Dfirst1(j,mm) = 1;
    end
    
    if (dtestfirst2(j) == d1)
        Dfirst2(j,mm) = 1;
    end
    
        if (dtestsecond1(j) == d2)
        Dsecond1(j,mm) = 1;
    end
    
    if (dtestsecond2(j) == d2)
        Dsecond2(j,mm) = 1;
    end
    
    end
end

mean(Dfirst1, 1)
mean(Dfirst2, 1)
mean(Dsecond1, 1)
mean(Dsecond2, 1)