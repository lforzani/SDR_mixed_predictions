function df = dF4newlad(W,FParameters)
sigma = FParameters.sigma;
sigmag = FParameters.sigmag;
nj = FParameters.n;
n = sum(nj);
p = cols(sigmag);

PHI = W;
A = FParameters.A;
ALPHA = PHI*A;
term1 = 2*inv(sigmag)*PHI*inv(PHI'*inv(sigmag)*PHI);
term2 = 2*sigmag*PHI*inv(PHI'*sigmag*PHI);
term3 = 2*sigmag*ALPHA*inv(ALPHA'*sigmag*ALPHA)*A';

f = nj/n;
a = zeros(rows(W), cols(W), length(f));
sigma_i = zeros(p);
for i=1:length(f)
    sigma_i(:,:) = sigma(i,:,:);
    a(:,:,i) = f(i)*sigma_i*ALPHA*inv(ALPHA'*sigma_i*ALPHA)*A';
end
  
df = term1 + term2 - term3 + 2*sum(a,3);
   



