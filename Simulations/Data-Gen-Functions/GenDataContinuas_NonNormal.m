function X = GenDataContinuas_NonNormal(fycent,Delta,A)

%Generaci?n de datos para el caso de predictores continuas
%r = 2; d=1, %dos predictores

p = size(Delta, 2);
n = size(fycent,1);
 
X = zeros(n,p);
for i = 1:n
  mu(i,:) = A*fycent(i,:)';
  a=(sqrtm(Delta));
  X(i,:) =mvtrnd(diag(ones(p,1)),6)*a + mu(i,:);
end


