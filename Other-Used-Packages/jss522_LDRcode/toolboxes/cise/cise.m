function [fv,beta,st] = cise(y,x,di,pw,method)
% This is the main function that achives coordinate-independent sparse estimation using
% either SIR or PFC for calculating M and N matrices. See references for details.
%    USAGE:
%  - outputs:
%    fv: the value of the objective function at the estimator
%    beta:  the estimator of the central subspace based on BIC criterion.
%    st: a vector with either 0 or 1 element to indicate which variable is
%    to be selected by CISE.

%  - inputs:
%    y: response vector.
%    x: predictors matrix.
%    di: the dimension of the central subspace
%    pw: the range of the penalty parameter of \lambda is from 0 to pw. 
%    method: model to calculate M and N matrices. So far, only 'PFC' or 'SIR' is accepted.
%
%
global p; % the number of predictor
global d; % the dimension of the central subspace
global n;
global M;
global N;

n = size(x,1);
p = size(x,2);
d= di;
nslice=6;

if strcmpi(method,'SIR')
    [M,N]=SIR(y,x,nslice);
elseif strcmpi(method,'PFC')
    [M,N] = PFC(y,x);
else
    error('Unknown method. Check syntax');
    exit;
end

nin = pw*100+1;

for i=1:nin
    f = biccis(i/100-0.01);
    stop(i) = f;
end

[f,ind] = min(stop);

lap = ind/100-0.01;

for k = 1:19
    f= biccis(lap-0.01+k/1000);
    dstop(k) = f;
end

[f,dind] = min(dstop);

dlap = lap-0.01+dind/1000;

[fv, beta, st] = biccis(dlap);
