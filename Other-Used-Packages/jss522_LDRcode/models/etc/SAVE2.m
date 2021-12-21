function [WX,W,auxmtx] = SAVE(Yaux,X,morph,dim,varargin)
% function [WX,W] = SAVE(Y,X,morph,u,varargin)
% 
% This function implements the sliced average variance estimation procedure
% for sufficient dimensionality reduction. See references for details.
%
% USAGE:
%   - outputs:
%     WX: projection of the predictors onto the dimension reduction subspace.
%     W: generating vectors of the dimension reduction subspace.
%     auxmtx: an auxiliary matrix with partial computations that can be
%     used to compute SIR and DR reductions.
%   - inputs:
%     Y: response vector.
%     X: predictors matrix.
%     morph: with value 'cont', specifies that the response Y is continuous 
%     (in which case it is a regression problem) while with value 'disc' it
%     specifies a discrete response (and a classification problem).
%     u: dimension for the reduced subspace.
%     varargin: optional arguments. They must be given as a 'OptionName', 'OptionValue' 
%     pair. Available options are limited to:
%     - 'nslices': to set the number of slices to be used to discretize continuous 
%       responses.
%     - 'auxmtx': to pass an auxiliary matrix computed with function
%     SETAUX, from which the SAVE reduction can be easily extracted. This is
%     useful when SAVE is used along with SIR and DR, as they share some
%     partial computations.

% =========================================================================

%----checking required arguments...........................................
if nargin < 4,
    error('Not enough input arguments. Type >>help ldr for details');
end
% .........................................................................

%----checking data consistency.............................................
if size(Yaux,1)~=size(X,1),
    error('The number of rows in Y must equate the number of rows in X');
end
if ~strcmpi(morph,'cont') && ~strcmpi(morph,'disc'),
    error('unknown type of response. Valid options are CONT or DISC...');
end
if ~ischar(dim) && ~isinZ(dim),
    error('Natural value expected to specify reduced subspace dimension.');
end
if ~isreal(Yaux),
    error('Response vector must be numeric');
end

%----reading optional input arguments and saving parameters................
parameters = read_input_nldr(varargin{:});

%----checking type of response and slicing if needed.......................
if strcmpi(morph,'disc'),
    Y = mapdata(Yaux);
    parameters.nslices = max(Y);
else % morph = 'cont'
    if parameters.nslices==0,
        warning('MATLAB:slices','for continuous responses, a number of slices should be given. Five slices will be used');
        parameters.nslices = 5;
    end
    Y = slices(Yaux,parameters.nslices);
end        

h = parameters.nslices;
[sigmas,means,sizes] = get_pars(X,Y,h);
delta = zeros(size(X,2));
for j=1:h,
    delta = delta + squeeze(sigmas(j,:,:));
end
% 
% 
% sigmag = get_cov(X);
[n,p] = size(X);
% [V,D] = eig(sigmag);
[V,D] = eig(delta);
if isempty(parameters.auxhandle),
    parameters.auxhandle = setaux(X,Y,parameters.nslices,p,V,D);
end
auxmtx = parameters.auxhandle;
[SAVEdir,SAVE] = getSAVE(dim,p,parameters.auxhandle);

%-----Write results------------------------
W = SAVE;
WX = (X - kron(ones(n,1),mean(X)))*W;
    
    
