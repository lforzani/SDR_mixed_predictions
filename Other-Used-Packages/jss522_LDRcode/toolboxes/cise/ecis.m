function [beta,st,k] = ecis(x)

global lambda;
global p;
global n;
global d;
global M;
global N;

if lambda ==0
    beta = x;
    st = ones(p,1);
    return;
end

d = size(x,2);

b0=x;
the = mywt(x);
%the = ones(p,1);
dis =1;
st = ones(p,1);
H0 = diag(the)*diag(1./sqrt(diag(b0*b0')));


[A,B] = mysqrt(N);
G = B*M*B;
k=0;

while dis > 1e-6

    k=k+1;
    [td,ty] = eig(G-0.5*lambda*B*H0*B);
    tp = length(H0);
    [a1,b1] = sort(diag(ty));
    if d==2
    b = B*[td(:,b1(tp)),td(:,b1(tp-1))];
    elseif d==1
        b = B*td(:,b1(tp));
    end
            
    dis = subspace(b0,b);
    st0 = st;
    [b,st,ix] = delp(b,st); 
    if sum(st) < sum(st0)    
        the(ix)=[];
        [A,B] = mysqrt(N);
        G = B*M*B;
    end;
    H0 = diag(the)*diag(1./sqrt(diag(b*b')));
    b0=b;
       if (size(H0,1)<=2) 
       break;
       end
end

te =0;

for i=1:p
    if st(i)==1
        te = te+1;
    beta(i,1:d) = b(te,1:d);
    elseif st(i)==0
        beta(i,1:d) = zeros(1,d);
    end
end
