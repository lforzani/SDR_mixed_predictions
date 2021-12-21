function  Vrcl = AsymptoticVariancebinarias(Tau0,Tauhat,fycent,q)
% RODRIGO: ver los argumentos que tenemos que pedir

% Matriz de varianza asintotica para el test de dimensi?n

n = size(fycent,1);
r = size(fycent,2);
m = (q+1)*q/2;
k = min(m,r);

%if (nargin == 6)
%    muH = zeros(q,1);
%end

c4 = zeros(q, r*q*(q-1)/2);
c6 = zeros(q*(q-1)/2,r*q);

Im = eye(m);
Iq = eye(q);
Ir = eye(r);
Iqr = eye(q*r);
Irq = Iqr;
Irqq = eye(r*q*(q-1)/2);
Iqq = eye(q*(q-1)/2);   

for i = 1:n
   faux = fycent(i,:); 
 
   e2 = kron(faux',Iq);
   e3 = kron(faux',Iqq);
   Fy = [Iq e2'  zeros(q,(r+1)*q*(q-1)/2); zeros(q*(q-1)/2,q + r*q) Iqq e3'];
   
   Derivada = derivbinaria2(Tau0,Tauhat,faux,q);
   
   auxV(i,:,:) = Fy'*Derivada*Fy;
end

mm = size(Fy,2);
au = mean(auxV,1);
invV = reshape(au(1,:,:) ,[mm,mm]);

V = inv(invV); 

M22 = [zeros(q*r,q) Iqr zeros(q*r, q*(q-1)/2+r*q*(q-1)/2)];
M25 = [zeros(r*q*(q-1)/2, q+q*r+q*(q-1)/2) Irqq];

M = [M22; M25];

Avarbtilde = M*V*M';

%W11 = [Ip; zeros(q,p); zeros(q*(q-1)/2,p)];
%W1 = kron(Ir,W11);
W2 = [Irq; zeros(r*q*(q-1)/2,r*q)];
%W2 = kron(Ir, W22);
W3 = [zeros(r*q,r*q*(q-1)/2); Irqq];
%W3 = kron(Ir, W33);

             W22=[Iq;zeros(q*(q-1)/2,q)];
             W33=[zeros(q,q*(q-1)/2);Iqq];
             
              W2=kron(Ir,W22);
             
             W3=kron(Ir,W33);
W = [W2 W3];
             
            
             

Vrcl = W*Avarbtilde*W';
   
end



