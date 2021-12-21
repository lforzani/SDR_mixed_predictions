function [red_estDlambda,T,redu] = rr4bin(H,Y,d,type,lambda)

% This function the minimal sufficient reduction for the regression 
% Y ~ H for binary data. 

% INPUTS:
% H: Binary predictor matrix, each row is an observation and each column a predictor variable.
% Y: response vector.
% d: dimension of the reduction 
% lambda: regularization parameter for variable selection.  Allowed values
% are 'auto' for authomatic variable selection, or a number indicanting the value of this parameter. 

% OUTPUTS:
% red_estDlambda: reduced-rank basis matrix for a dimension reduction subspace 
% for binary data with variable selection. 
% T: sufficient statistics. 
% redu: reduced-rank basis matrix for a dimension reduction subspace 
% for binary data without variable selection.

%%
if nargin < 4,
    error('not enough input parameters');
end
if nargin < 5,
    lambda = 0.0;
end

[n,q] = size(H);
if length(Y)~=n,
    error('number of rows in H and Y must be the same')
end

if isstr(lambda),
    [T,~,C] = EM4mixture_binaria(H,Y,type);
    [p] = size(C,1);
    [A,B] = initAB(C,d);
    opts = initOpts4bin(B',C',A');
    [XX,YY] = getXY(B',C');
    lambda_max = 2*getLambdaMax4bin(YY,XX,opts);
    
    % grid
    nlambda = 50;
    z1 = logspace(0,-2,nlambda);
    lambdas = z1*lambda_max;

    % cross validation
    kfold = 8;
    errores = ones(kfold,nlambda);
    fold = randi(kfold,n,1);
    for k=1:kfold
        idx = find(fold==k);
        datos.Xtrain = H; 
        datos.Xtrain(idx,:) = [];
        datos.Xtest = H(idx,:);
        datos.Ytrain = Y; 
        datos.Ytrain(idx) = [];
        datos.Ytest = Y(idx);
        [~,~,redu] = EM4mixture_binaria(datos.Xtrain,datos.Ytrain,type);
        [~,errores(k,:)] = getLambda4bin(redu,d,datos,lambdas);
    end
    [~,kopt] = min(mean(errores));
    lambda = lambdas(kopt);
end

if lambda > 0,
    [T,~,redu,proj,fycent] = EM4mixture_binaria(H,Y,type);
    [red_estDlambda] = getABpen4bin_simple(redu,d,lambda);
else
    [T,~,redu,proj,fycent] = EM4mixture_binaria(H,Y,type);
    [red_estDlambda] = initAB(redu,d);
end
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
