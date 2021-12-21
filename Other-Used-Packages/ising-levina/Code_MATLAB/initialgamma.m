% Initialize gamma using penalized logistic regression
% 
% Call    : a = initialgamma(yk,yj, x,lambda)
%
% Arguments
% yk      : k-th column of y matrix
% yj      : j-th column of y matrix
% x       : n x p covariate matrix
% lambda  : penalty parameter
%
% Values
% a       : returned scalar initial value

function[a] = initialgamma(yk,yj, x,lambda)

%% input yk yj are n-dim 0-1 vector, x is n-dim cts vector.
y = [yk ; yj];
z = [x.*yj ; x.*yk];

a = initialbeta(y,z,lambda);

end