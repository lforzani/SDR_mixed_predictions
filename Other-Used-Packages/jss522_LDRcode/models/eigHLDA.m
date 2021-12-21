function [WX,W] = eigHLDA(Yaux,X,u,morph)
%----checking type of response ......................
if strcmpi(morph,'disc'),
    Y = grp2idx(Yaux);
    parameters.nslices = max(Y);
else % morph = 'cont'
    Y = Yaux;
    parameters.nslices = length(Y);
end        

%--- get sample statistics ................................................
data_parameters = setdatapars(Y,X,parameters.nslices);
SIGMA = data_parameters.sigmag;
SIGMAS = data_parameters.sigma;
n=data_parameters.n;

Shlda = -sum(n)*logmat(SIGMA);
% sigma_i = zeros(size(SIGMA));
for i=1:size(SIGMAS,1),
    sigma_i = squeeze(SIGMAS(i,:,:));
    Shlda= Shlda + n(i)*logmat(sigma_i);
end
Shlda = -Shlda/2;

[V,D]=firsteigs(Shlda,u);
W=orth(V);
WX=X*W;

% function R = logmat(M)
% %logarithm of a matrix
% [V D]= eig(M);
% d = diag(D);
% logd = log(d);
% D = diag(logd);
% R = V * D * V';
