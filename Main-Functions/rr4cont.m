function [red_estDlambda,T,redu] = rr4cont(X,Y,d,type,lambda)

% This function the minimal sufficient reduction for the regression 
% Y ~ X for continuous data. 

% INPUTS:
% X: Continuos predictor matrix, each row is an observation and each column a predictor variable.
% Y: response vector.
% d: dimension of the reduction 
% type: string or number used to set the type of the response. Allowed values are 'disc' for discrete 
% responses, or a number indicating the degree of the polinomial of Y. 
% lambda: regularization parameter for variable selection.  Allowed values
% are 'auto' for authomatic variable selection, or a number indicanting the value of this parameter. 

% OUTPUTS:
% red_estDlambda: reduced-rank basis matrix for a dimension reduction subspace 
% for continuous data with variable selection. 
% T: sufficient statistics. 
% redu: reduced-rank basis matrix for a dimension reduction subspace 
% for continuous data without variable selection.

%% 
if nargin < 4,
    error('not enough input parameters');
end
if nargin < 5,
    lambda = 0.0;
end

[n,p] = size(X);
if length(Y)~=n,
    error('number of rows in X and Y must be the same')
end

if isstr(lambda),
    [~,C] = EM4mixture_continua(X,Y,type);
    [p] = size(C,1);
    [A,B] = initAB(C,d);
    opts = initOpts(A');
    [XX,YY] = getXY(B',C');
    lambda_max = getLambdaMax4cont(YY,XX,opts);
    
    % grid
    nlambda = 100;
    z = logspace(0,-5,nlambda);
    lambdas = z*lambda_max/50;

    % cross validation 
    kfold = 5;
    errores = ones(kfold,nlambda);
    fold = randi(kfold,n,1);
    for k=1:kfold
        idx = find(fold==k);
        datos.Xtrain = X; 
        datos.Xtrain(idx,:) = [];
        datos.Xtest = X(idx,:);
        datos.Ytrain = Y; 
        datos.Ytrain(idx) = [];
        datos.Ytest = Y(idx);
        [~,redu] = EM4mixture_continua(datos.Xtrain,datos.Ytrain,type);
        [~,errores(k,:)] = getLambda4cont(redu,d,datos,lambdas);
    end
    [~,kopt] = min(mean(errores));
    lambda = lambdas(kopt);
end
if lambda > 0,
    [T,redu,proj,Ahat,Deltahat,fycent] = EM4mixture_continua(X,Y,type);
    [red_estDlambda] = getABpen_vsimple(redu,d,lambda);
else
    [T,redu,proj,Ahat,Deltahat,fycent] = EM4mixture_continua(X,Y,type);
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

function [X,Y] = getXY(A,C)
    r = size(C,2);
    kronn = kron(eye(r),A);
    X = kronn; 
    Y = C(:);
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

     
