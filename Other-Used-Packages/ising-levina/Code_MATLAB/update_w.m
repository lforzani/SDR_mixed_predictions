% Function to update the weight matrix when probability matrix is updated
%
% Call  : wnew = update_w(w,prob,j,k)
%
% Arguments
% w     : n x q old weight matrix
% prob  : n x q updated probability matrix
% j     : 2nd co-ordinate of the updated gamma parameter
% k     : 3rd co-ordinate of the updated gamma parameter
%
% Values
% wnew  : n x q updated weight matrix


function[wnew] = update_w(w, prob, j, k)

wnew = w;
if(j ~= k)
    wnew(:,j) = prob(:,j).*(1-prob(:,j));
end
wnew(:,k) = prob(:,k).*(1-prob(:,k));
end
