
%load('Y620Lasso_dummy.mat');
q = 620;
p = 4;
L = length(lambda);
for l = 1:L
    l
for i = 1:p
    name = sprintf('DataAnalysis/Y%d_ENMin_X%d_lambda%d.csv',q,i,l);
    tmp = reshape(gamma_min(i+1, :,:,l),q,q);
    csvwrite(name, tmp);
end
end





