function f = F4hlda(W,FParameters)
% 
% W: projection matrix onto the reduced subspace
% FParameters: structure with needed statistics

sigmas = FParameters.sigma;
SIGMA = FParameters.sigmag;
nj = FParameters.n;
n = sum(nj);
p = cols(SIGMA);
% ---define some convenience variables
h = nj/n;

a = zeros(length(h),1);
sigma_i = zeros(p);
for i=1:length(h),
    sigma_i(:,:) = sigmas(i,:,:);
    a(i) = logdet(W'*sigma_i*W);
end

f = n/2*(logdet(W'*inv(SIGMA)*W) + h*a);