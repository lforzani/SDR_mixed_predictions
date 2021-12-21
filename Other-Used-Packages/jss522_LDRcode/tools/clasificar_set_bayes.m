function [error,labels] = clasificar_set(testset,betas,modelos,nclases,priors)
Ytest = testset(:,1); Xtest = testset(:,2:end);
labels = zeros(size(Ytest));
for n=1:length(Ytest),
    labels(n) = clasificar_one_bayes(Xtest(n,:),betas,modelos,nclases,priors);
end
error = length(find(labels~=Ytest))/length(Ytest);


