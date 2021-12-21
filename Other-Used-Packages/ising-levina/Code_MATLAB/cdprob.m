% Function to get the conditional probability matrix
%
% Call  : prob = cdprob(x,y,par)
%
% Arguments
% x     : n x p covariate matrix
% y     : n x q binary response matrix
% par   : (p+1) x q x q dimensional parameter array
%
% Values
% prob  : n x q dimensional matrix of conditional probabilities of 
%         Pr(y_j = 1 | x, y_(-j))


function[prob] = cdprob(x,y,par)
x1 = [ones(size(x,1),1) x];
q = size(y,2);

for j = 1:q
    A = x1*par(:,:,j);
    z = y;
    z(:,j) = 1;
    for k = 1:size(x,1)
        prob(k,j) = 1/(1+exp(-A(k,:)*z(k,:)'));
    end
%     ita_t = A*z';
%     %ita(:,j) = diag(ita_t);
%     prob(:,j) = 1./(1+exp(-diag(ita_t)));
end

end