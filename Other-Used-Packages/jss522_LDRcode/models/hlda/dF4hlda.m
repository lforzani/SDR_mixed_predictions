function df = dF4hlda(W,FParameters)

sigmas = FParameters.sigma;
SIGMA = FParameters.sigmag;
nj = FParameters.n;
n = sum(nj);
p = cols(SIGMA);
% ---define some convenience variables
h = nj/n;

sigma_i = zeros(p);
df = zeros(size(W));
for i=1:length(h),
    sigma_i(:,:) = sigmas(i,:,:);
    df = df + h(i)*sigma_i*W*inv(W'*sigma_i*W);
end

df = n*(df + inv(SIGMA)*W*inv(W'*inv(SIGMA)*W));