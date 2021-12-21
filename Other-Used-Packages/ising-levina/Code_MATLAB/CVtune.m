% function to choose the tuning parameter by Cross Validation

function[lambda_opt likeli] = CVtune(y, x, lambda, nfold, option)
n = size(y,1);
L = length(lambda);

ind = crossvalind('Kfold',n,nfold);
likeli = zeros(L,nfold);
for k = 1:nfold
    xtrain = x(find(ind==k),:);
    ytrain = y(find(ind==k),:);
    xtest = x;
    ytest = y;
    xtest(find(ind==k),:)=[];
    ytest(find(ind==k),:)=[];
    tmp = zeros(L,1);
    gamma_mat = SolPath(ytrain, xtrain, lambda, option);
    if strcmp(option, 'separate')
        gamma_mat = select_gamma(gamma_mat);
    end
    for l = 1:L
    tmp(l) = loss(ytest,xtest,gamma_mat(:,:,:,l), 0);
    end
    likeli(:,k) = tmp;
end
[~, lambda_opt] = min(mean(likeli,2));
end

