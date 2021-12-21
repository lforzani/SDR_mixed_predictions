function [WTaida,Taida] = aida(Y,X,u,morph)
if nargin<4,
    morph='disc';
end

%the size of the data space  
[Nt p] = size(X);

%the number of classes
if strcmpi(morph,'cont'),
    Y = slices(Y,5);
end
Y = grp2idx(Y);
h = max(Y);

data_parameters = setdatapars(Y,X,h);
prob_i = data_parameters.n/sum(data_parameters.n);
Sigma_i = data_parameters.sigma;
Sigma = data_parameters.sigmag;
Sigma_w = zeros(p);
for i=1:h,
    Sigma_w = Sigma_w + prob_i(i)*squeeze(Sigma_i(i,:,:));
end
% Sigma_b = Sigma - Sigma_w;

% TRANSFORM0
%sphering
rootSw = sqrtmat(Sigma_w);   %transforms the variables so Sigma_w = I;
W = inv(rootSw);             %whitening matrix so that Sigma_w = I;
% Zi=zeros(p,p,h);
% for i = 1:h
%     Zi(:,:,i) = W * squeeze(Sigma_i(:,:,i)) * W;  %transformed si's
% end 
S_T = W*Sigma*W;
Saida = logmat(S_T);
for i=1:h,
    Saida = Saida - prob_i(i)*logmat(W * squeeze(Sigma_i(i,:,:)) * W);
end

[V,D] = firsteigs(Saida,u);
%Taida = (V/norm(V));
Taida = orth(V);
WTaida = W*Taida;

function R = logmat(M)
%logarithm of a matrix

[V D]= eig(M);
d = diag(D);
logd = log(d);
D = diag(logd);
R = V * D * V';

function R = sqrtmat(M)
%square root of a matrix
%note that root is a symmetric matrix

[V D]= eig(M);
R = V * D.^0.5 * V'; 
