function f = F4newlad(W,FParameters)
sigma = FParameters.sigma;
sigmag = FParameters.sigmag;
nj = FParameters.n;
n = sum(nj);
p = cols(sigmag);
PHI = W;
A = FParameters.A;

% ---define some convenience variables
term1 = logdet(PHI'*inv(sigmag)*PHI);
term2 = logdet(PHI'*sigmag*PHI);
ALPHA = PHI*A; 
term3 = logdet(ALPHA'*sigmag*ALPHA);

h = nj/n;
a = zeros(length(h),1);
sigma_i = zeros(p);
for i=1:length(h),
    sigma_i(:,:) = sigma(i,:,:);
    a(i) = logdet(ALPHA'*sigma_i*ALPHA);
end

% ---Likelihood function for LAD model
f = term1 + term2 - term3 + h*a;

