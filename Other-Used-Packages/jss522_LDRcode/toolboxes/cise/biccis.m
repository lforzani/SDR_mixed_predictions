function [f,beta,st] = biccis(lam)

global lambda;
global n;
global d;
global p;
global M;
global N;  

M0= M;
N0= N;

if lam < 0
    f = inf;
    return
end

lambda = lam;

[A,B] = mysqrt(N);
G = B*M*B;
[a2,b2]= eig(G);
[a1,b1] = sort(diag(b2));

if d == 2
x0 = B*[a2(:,b1(p)),a2(:,b1(p-1))];

elseif d==1
    x0 = B*a2(:,b1(p));
end

[beta,st] = ecis(x0);

M = M0;
N = N0;



%f =  -trace(beta'*M*beta)+ 2*d*(sum(st)-d)/n;
f = -trace(beta'*M*beta) + d*(sum(st)-d)*log(n)/n;
%f = p*log(trace(Sigma0)-trace(beta'*Sigma0*beta)) + 2*d*sum(st)*2/n;

%RIC
%pe = d*(sum(st)-d);
%f = -(n-pe)*trace(beta'*M*beta) + (n-pe)*log(n/(n-pe)) + pe*(log(n)-1) + 4/(n-pe-2);
