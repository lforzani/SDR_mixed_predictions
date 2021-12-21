

n = 40;
x1 =  normrnd(0, 1, n, 1);
x2 = 0.5*x1+normrnd(0,5,n,1);
x3 = mvnrnd(zeros(20-2,1), 0.5*eye(20-2), n);
logit = x2+x1+normrnd(0, 0.05, n, 1);
y = binornd(1, 1./(1+exp(-logit)));
x =  [x1 x2 x3];
lambda = 0.5.^(0:1:5);

options = glmnetSet;
options.lambda = lambda;
options.standardize = 0;

freq = zeros(size(x,2), length(lambda));
for k = 1:1000
    k
    ind  = randsample(n, n/2);
fit = glmnet(x(ind,:),y(ind, :)+1,'binomial',options);
freq = freq+(fit.beta~=0);
end
freq_max = max(freq,[],2)
