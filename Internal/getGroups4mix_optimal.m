function [g,w] = getGroups4mix_optimal(p,q,d)
% estrategia: hacemos los grupos para las binarias como antes y luego
% pegamos adelante la parte para las continuas...
%
%
% grupos para la parte binaria

pp = p+q*(q+1)/2;

a = magic(q);
at = tril(a)-diag(diag(a));
idx=find(at);
at(idx)=1:length(idx);
for dd=1:(q-1),
    at = at + diag(diag(at,-dd),dd);
end
g = at(at>0);
g = [1:q,g'+q];
ng = length(g);
 
gbin = g;
aux = pp-p;
for dd=2:d,
    gbin = [gbin; g+(dd-1)*aux];
end
% for dd=2:d,
%     gbin = [gbin; g+(dd-1)*ng];
% end

total = q*q*d;
starts1 = 1:d:((q-1)*d+1);
ends1 = d:d:q*d;
starts2 = (q*d+1):(q-1)*d:(total-(q-1)*d+1);
ends2 = (q*d+(q-1)*d): (q-1)*d : total;

% juntamos las dos partes y sumamos p, ya que hay p grupos asociados a las
% variables continuas
starts_bin = [starts1,starts2] + d*p;
ends_bin =  [ends1,ends2] + d*p;
% modificos los n?meros de las variables en los grupos, para incluir las
% continuas

gmix = [zeros(d,p),gbin];
paux = size(gbin,2)+p;

for dd=1:d,
    gmix(dd,1:p) = ((dd-1)*pp + 1):((dd-1)*pp + p);
    gmix(dd,(p+1):end) = gbin(dd,:)+dd*p;
end

% for dd=1:d,
%     gmix(dd,1:p) = ((dd-1)*pp + 1):(dd*p);
%     gmix(dd,(p+1):end) = gbin(dd,:)+dd*p;
% end

g = gmix(:); g=g';     

w = ones(3,p+2*q);
w(1,:)=[1:d:(p-1)*d+1,starts_bin];
w(2,:)=[d:d:p*d,ends_bin];
w(3,:)= sqrt(w(2,:)-w(1,:) + 1);
end
