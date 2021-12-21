function [g,w] = getGroups4bin(p)
a = magic(p);
at = tril(a)-diag(diag(a));
idx=find(at);
at(idx)=1:length(idx);
for d=1:(p-1),
    at = at + diag(diag(at,-d),d);
end
g = at(at>0);
starts = 1:(p-1):(length(g)-(p-1)+1);
ends = (p-1):(p-1):length(g);
w = ones(3,2*p);
w(1:2,:)=[1:p p+starts; 1:p, p+ends];
w(3,:) = [ones(1,p) sqrt(p-1)*ones(1,p)];
 g = [1:p,g'+p];