function [WX,W] = eigEPFC(Yaux,X,u,morph,FF)
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
%the size of the data space  
[Nt p] = size(X);

if strcmpi(morph,'cont')
    SIGMAfit = get_fitted_cov(Y,X,FF);
else
    SIGMAfit = get_average_cov(X,data_parameters);
end
SIGMA = data_parameters.sigmag;
SIGMAres = SIGMA - SIGMAfit;


% S = -.5*(logm(inv(SIGMA))+logm(SIGMAres));
S = .5*(logmat(SIGMA)-logmat(SIGMAres));
[W,D] = firsteigs(S,u);
W=orth(W);
WX=X*W;


function R = logmat(M)
%logarithm of a matrix

[V D]= eig(M);
d = diag(D);
logd = log(d);
D = diag(logd);
R = V * D * V';
