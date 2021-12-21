% Function to compute the soft-thresholding of a given vector with respect
% to a thresholding vector of same length
% 
% Call : s = lasso(y,r)
%
% Arguments
% y    : p x 1 vector to e softly shrunk
% r    : p x 1 vector of threshold parameter
%
% Values
% s    : p x 1 vector solution

function[s] = lasso(y,r)

p = size(y,1);
for i = 1:p
if       y(i)>0 && r(i)<abs(y(i))
                           t(i) = y(i)-r(i);
elseif   y(i)<0 && r(i)<abs(y(i))
                           t(i) = y(i)+r(i);
else
                           t(i) = 0;
end
end
s = t';
end