function [A] = rr4mix_optimal(X,H,Y,d,type,lambda,thr,sparsepar)

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
% A: reduced-rank basis matrix for a dimension reduction subspace 
% for mixed (binary and continuous) data.

%%

if nargin < 5,
    error('not enough input parameters');
end
if nargin < 6,
    lambda = 0.0;
end
if nargin < 7
    thr = 1e-3;
end

if nargin < 8
    sparsepar = 0.6;
end


[n,q] = size(H);
[m,p] = size(X);
if n~=m,
    error('number of rows in H and Y must be the same');
end
if length(Y)~=n,
    error('number of rows in H and Y must be the same')
end

nlambda = 50;
[C] = EM4mixture_continua_binaria(X,H,Y,type);
[pp] = size(C,1);

if isstr(lambda),
    ppp = (p+q*(q+1)/2);

    [A,B] = initAB(C,d);
    opts = initOpts(B',C',A',pp,p,q);
    [XX,YY] = getXY(B',C');
    opts.d = d;
    lambda_max = getLambdaMax4cont_optimal(YY,XX,opts);
    
    % grid
    z1 = logspace(-1,-3,nlambda);
    lambdas = z1*lambda_max;

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
        [redu] = EM4mixture_continua_binaria(datos.Xtrain(:,1:p),datos.Xtrain(:,(p+1):(p+q)),datos.Ytrain,type);
        [errores(k,:)] = getLambda4optimal_v2(redu,d,p,q,datos,lambdas);
    end
    [~,kopt] = min(mean(errores));
    lambda = lambdas(kopt);
end

%%
[A,B] = initAB(C,d);

if lambda > 0,
    A = updateB(B',C',A',pp,p,q,lambda,sparsepar); A = A';
    A=orth(A);
    for j=1:pp,
        if norm(A(j,:)) < thr,
            A(j,:)=0;
        end
    end
end
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

function B = updateB(A,C,B0,pp,p,q,lambda,sparsepar)
% esta funcion acutaliza B utilizando el metodo iterativo
% de minimizacion cuadratica presentado en la pagina 47
% de la tesis de sabri
[~,r] = size(C);
d = size(A,2); 

if size(B0,2) > 1,
    B0 = B0(:);
end

if lambda == 0,
    kronn = kron(eye(r),A);
    aux = inv(kronn'*invV*kronn)*kronn'*invV*C(:);
else
    opts=[];
    opts.init=B0;        % starting point
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
    
    z = [sparsepar*lambda, (1-sparsepar)*lambda];
    %z = [0 lambda]; 

    [G,w] = getGroups4mix_optimal(p,q,d);
    opts.G = G;
    opts.ind = w;
    opts.x0 = B0;

    kronn = kron(eye(r),A);
    X = kronn; 
    Y = C(:);
    [aux, funVal3, ValueL3]= overlapping_LeastR(X, Y, z, opts);

end    

B = reshape(aux,[d,r]);
end


function out = sqrm(A)
[V,d] = eig(A);
out = V*diag(sqrt(d))*V';
end


function opts = initOpts(A,C,B0,pp,p,q)
[~,r] = size(C);
d = size(A,2); 

if size(B0,2) > 1,
    B0 = B0(:);
end

    opts=[];
    opts.init=B0;        % starting point
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
    
    [G,w] = getGroups4mix_optimal(p,q,d);
    opts.G = G;
    opts.ind = w;
    opts.x0 = B0;
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

function [X,Y] = getXY(A,C)
    r = size(C,2);
    kronn = kron(eye(r),A);
    X = kronn; 
    Y = C(:);
end


