function [red_est_suboptimal1,red_est_suboptimal2,betahat] = rr4mix_suboptimal(X,H,Y,d1,d2,type,lambdaX,lambdaH)

% This function the minimal sufficient reduction for the regression 
% Y ~ (X,H) for binary and continuous data. 

% INPUTS:
% X: Continuos predictor matrix, each row is an observation and each column a predictor variable.
% H: Binary predictor matrix, each row is an observation and each column a predictor variable.
% Y: response vector.
% d: dimension of the reduction 
% type: string or number used to set the type of the response. Allowed values are 'disc' for discrete 
% responses, or a number indicating the degree of the polinomial of Y. 
% lambda: regularization parameter for variable selection.  Allowed values
% are 'auto' for authomatic variable selection, or a number indicanting the value of this parameter. 
% sparsepar,thr: hyperparameters for parameter initialization. 

% OUTPUTS:
% red_est_suboptimal1: reduced-rank basis matrix for a dimension reduction subspace 
% for the continuous part of the reduction.
% red_est_suboptimal2: reduced-rank basis matrix for a dimension reduction subspace 
% for the binary part of the reduction. 
% betahat: estimated parameter beta (see the paper), neccessary to compute the projection 
%of continuous part of the suboptimal reduction.


% A: reduced-rank basis matrix for a dimension reduction subspace 
% for mixed (binary and continuous) data.

%% 
if nargin < 6,
    error('not enough input parameters');
end
if nargin < 7,
    lambdaX = 0.0;
end
if nargin < 8,
    lambdaH=0.0;
end

[n,q] = size(H);
[m,p] = size(X);
if n~=m,
    error('number of rows in H and Y must be the same');
end
if length(Y)~=n,
    error('number of rows in H and Y must be the same')
end


[~,redu_suboptimal, fycent,T_optimal,T_suboptimal,Tauhatraro0, Tauhatraro, Ahat, betahat, Deltahat] = EM4mixture_continua_binaria(X,H,Y,type);
[pp] = size(redu_suboptimal,1);
C1 = redu_suboptimal(1:(p+q),:);
C2 = redu_suboptimal(1+(p+q):pp,:);

if isstr(lambdaX),
    
    % (X U H)
    [A1,B1] = initAB(C1,d1);
    opts1 = initOpts(A1');
    [XX,YY] = getXY(B1',C1');
    lambdaX_max = getLambdaMax4cont(YY,XX,opts1);
    
    % grid
    nlambda = 50;
    z1 = logspace(-1,-3,nlambda);
    lambdas = z1*lambdaX_max;

    % cross validation 
    kfold = 5;
    errores = ones(kfold,nlambda);
    fold = randi(kfold,n,1);
    for k=1:kfold
        idx = find(fold==k);
        datos.Xtrain = [X,H]; 
        datos.Xtrain(idx,:) = [];
        datos.Xtest = [X(idx,:),H(idx,:)];
        datos.Ytrain = Y; 
        datos.Ytrain(idx) = [];
        datos.Ytest = Y(idx);
        [~,redu] = EM4mixture_continua_binaria(datos.Xtrain(:,1:p),datos.Xtrain(:,(p+1):(p+q)),datos.Ytrain,type);
        [~,errores(k,:)] = getLambda4cont(redu(1:(p+q),:),d1,getCont(datos,p+q),lambdas);
    end
    [~,kopt] = min(mean(errores));
    lambdaX = lambdas(kopt);
    
    
    %% para la porcion HH
    [qq] = size(C2,1);
    [A2,B2] = initAB(C2,d2);
    opts2 = initOpts4bin(B2',C2',A2');
    [XX,YY] = getXY(B2',C2');
    lambdaH_max = getLambdaMax4bin(YY,XX,opts2);
    
    % grid
    z1 = logspace(1,-3,nlambda);
    lambdas = z1*lambdaH_max;

    % cross validation
    kfold = 5;
    errores = ones(kfold,nlambda);
    fold = randi(kfold,n,1);
    for k=1:kfold
        idx = find(fold==k);
        datos.Xtrain = [X,H]; 
        datos.Xtrain(idx,:) = [];
        datos.Xtest = [X(idx,:),H(idx,:)];
        datos.Ytrain = Y; 
        datos.Ytrain(idx) = [];
        datos.Ytest = Y(idx);
        [~,redu] = EM4mixture_continua_binaria(datos.Xtrain(:,1:p),datos.Xtrain(:,(p+1):(p+q)),datos.Ytrain,type);
        [~,errores(k,:)] = getLambda4bin(redu(p+q+1:end,:),d2,getBin(datos,p),lambdas);
    end
    [~,kopt] = min(mean(errores));
    lambdaH = lambdas(kopt);
end    
    
    [A1,~] = getABpen_vsimple(C1,d1,lambdaX);
    [A2,~] = getABpen4bin_simple(C2,d2,lambdaH);

    red_est_suboptimal1 = A1;
    red_est_suboptimal2 = A2;

end


function opts = initOpts(B0)
    [d,r] = size(B0);
    opts=[];
    opts.maxIter=1000;   % maximum number of iterations
    opts.nFlag=0;       % without normalization
    opts.rFlag=0;       % the input parameter 'rho' is a ratio in (0, 1)
    opts.rsL2=0.0;     % the squared two norm term
    opts.mFlag=1;       % smooth reformulation 
    opts.lFlag=1;       % adaptive line search
    opts.tFlag=5; 
    opts.ind = [0, d*(1:r)];
    opts.x0 = B0(:);
end

function [Aini,Bini] = initAB(C,d)
% esta funcion calcula A y B iniciales utilizando el estimador
% Dub-d

[U,D,V] = svd(C);

Ud = U(:,1:d);
Vd = V(:,1:d); 
Dd = D(1:d,1:d);

Aini = Ud*Dd;
Bini = Vd';
end

     

function opts = initOpts4bin(A,C,B0)
    [p,d] = size(C);
    d = size(A,2);
    opts=[];
    opts.maxIter=5000;   % maximum number of iterations
    opts.maxIter2 = 1000;
    opts.nFlag=0;       % without normalization
    opts.rFlag=0;       % the input parameter 'rho' is a ratio in (0, 1)
    % the squared two norm term
    opts.mFlag=0;       % smooth reformulation 
    opts.lFlag=0;       % adaptive line search
    opts.tFlag=3;
    opts.tol = 1e-8;
    opts.tol2 = 1e-8;
    opts.flag2 = 2;
    [G,w] = getGroups4bin(p,d);
    opts.G = G;
    opts.ind = w;
    opts.x0 = B0(:);
end

function [X,Y] = getXY(A,C)
    r = size(C,2);
    kronn = kron(eye(r),A);
    X = kronn; 
    Y = C(:);
end


function T = getT(H)
[n,q] = size(H);
k = q*(q+1)/2;
T= zeros(n,k);
    for j=1:n
        auxt = H(j,:)'*H(j,:);
    T(j,:) = [diag(auxt);auxt(find(tril(ones(q,q),-1)))];
    end
end