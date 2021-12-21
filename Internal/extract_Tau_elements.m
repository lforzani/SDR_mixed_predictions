function [Lqtau, Jqtau] = extract_Tau_elements(tau)

k = size(tau,1);
q = (-1+ sqrt(1+8*k))/2;
%k2 = q*(q-1)/2;
r = size(tau,2);

pos = ones(1,q);
algo = 0;
for i = 1:(q-1)
    algo = algo + (i-1);
    pos(i+1) = (i*q - algo) + 1;
end 

Lqtau = tau(pos,:);
aux = tau;
aux(pos,:) = [];
Jqtau = aux;


