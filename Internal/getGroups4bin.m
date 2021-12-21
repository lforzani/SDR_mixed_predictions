function [g,w] = getGroups4bin1D(p,d)
a = magic(p);
at = tril(a)-diag(diag(a));
idx=find(at);
at(idx)=1:length(idx);
for dd=1:(p-1),
    at = at + diag(diag(at,-dd),dd);
end
g = at(at>0);
g = [1:p,g'+p];
%ng = length(g);
ng = p*(p+1)/2; 
gd = g;
for dd=2:d,
    gd = [gd; g+(dd-1)*ng];
end
total = p*p*d;
starts1 = 1:d:((p-1)*d+1);
ends1 = d:d:p*d;
starts2 = (p*d+1):(p-1)*d:(total-(p-1)*d+1);
ends2 = (p*d+(p-1)*d): (p-1)*d : total;

w = ones(3,2*p);
w(1,:)=[starts1,starts2];
w(2,:)=[ends1,ends2];
     w(3,:)= sqrt(w(2,:)-w(1,:) + 1);
g = gd(:); g=g';     
