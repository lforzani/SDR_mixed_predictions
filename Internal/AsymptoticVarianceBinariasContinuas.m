function  Vrcl = AsymptoticVarianceBinariasContinuas(Tau0,Deltahat,Ahat,Tauhat,betahat,fycent,muX,muH);

% Matriz de varianza asintotica para el test de dimensi?n


n = size(fycent,1);
p = size(Ahat,1);
q = size(betahat,2);
r = size(fycent,2);
m = p+(q+1)*q/2;
k = min(m,r);

if (nargin == 6)
    muX = zeros(p,1);
    
    muH = zeros(q,1);
end

c1 = zeros(p,q*r);
c2 = zeros(p,r*q*(q-1)/2);
c3 = zeros(q, r*p);
c4 = zeros(q, r*q*(q-1)/2);
c5 = zeros(q*(q-1)/2, r*p);
c6 = zeros(q*(q-1)/2,r*q);

Ip = eye(p);
Iq = eye(q);
Ir = eye(r);
Ipq = eye(p*q);
Ipr= eye(p*r);
Iqr =eye(q*r);
Irqq = eye(r*q*(q-1)/2);
Iqq = eye(q*(q-1)/2);   
Ipp = eye(p*(p+1)/2);

for i = 1:n
   faux = fycent(i,:)'; 
   
   e1 = kron(faux',Ip);
   e2 = kron(faux',Iq);
   e3 = kron(faux',Iqq);
   F1 = [Ip e1 zeros(p,q*(1+p+r+(1+r)*(q-1)/2)+ p*(p+1)/2)];
   F2 = [zeros(q,p*(1+r)) Iq e2 zeros(q,q*(p+(1+r)*(q-1)/2)+ p*(p+1)/2)];
   F3 = [zeros(p*(p+1)/2,(p+q)*(1+r)) Ipp zeros( p*(p+1)/2,q*(p+(1+r)*(q-1)/2))];
   F4 = [zeros(p*q,(p+q)*(1+r)+ p*(p+1)/2) Ipq zeros( p*q,q*(1+r)*(q-1)/2)];
   F5 = [zeros(q*(q-1)/2,(p+q)*(1+r)+ p*(q+(p+1)/2)) Iqq e3];
   Fy = [ F1; F2 ; F3 ; F4 ; F5] ;
   
   Derivada = derivadacompleta(muX,muH,Tau0,Deltahat,Ahat,betahat,Tauhat,faux);
    
   auxV(i,:,:) = Fy'*Derivada*Fy;
end

mm = size(Fy,2);
au = mean(auxV,1);
invV = reshape(au(1,:,:) ,[mm,mm]); % + 1*eye(mm);%ojo aca

%Aqui probamos con la diagonal por que es muy lento. En realidad va
%V=inv(invV)
%DinvV=diag(diag(invV));
%V = inv(DinvV); 

V=inv(invV);

M21 =[zeros(p*r,p) Ipr  zeros(p*r, q+q*r+p*(p+1)/2+p*q+q*(q-1)/2+r*q*(q-1)/2)];
M22 =[zeros(q*r,p+p*r+q) Iqr zeros(q*r, p*(p+1)/2+p*q+q*(q-1)/2+r*q*(q-1)/2)];
M25 =[zeros(r*q*(q-1)/2, p+p*r+q+q*r+p*(p+1)/2+p*q+q*(q-1)/2) Irqq];

M = [M21; M22; M25];

Avarbtilde = M*V*M';


W11=[ Ip; zeros(q,p); zeros(q*(q-1)/2,p)];
W1 = kron(Ir,W11);
W22=[ zeros(p,q); Iq; zeros(q*(q-1)/2,q)];
W2 = kron(Ir, W22);
W33=[ zeros(p,q*(q-1)/2); zeros(q,q*(q-1)/2); Iqq];
W3 = kron(Ir, W33);

W=[ W1 W2 W3];

Vrcl = W*Avarbtilde*W';
   
end



