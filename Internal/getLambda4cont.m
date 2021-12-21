function [lambda_opt,performance] = getLambda4cont(C,d,datos,nlambda)
if nargin < 4,
    nlambda = 20;
end

[p] = size(C,1);
[A,B] = initAB(C,d);
opts = initOpts(A');
[X,Y] = getXY(B',C');
if length(nlambda)==1,
    lambda_max = getLambdaMax4cont(Y,X,opts);
    z = logspace(1,-3,nlambda);
    lambdas = z*lambda_max;
else
    lambdas=nlambda;
end

Xtrain = datos.Xtrain;
Ytrain = datos.Ytrain;
Xtest = datos.Xtest;
Ytest = datos.Ytest;


[vecAt, funVal3, ValueL3]= glLeastR(X, Y, lambdas(1) , opts);
 A = reshape(vecAt,d,p); A = A';
 nlambda = length(lambdas);
performance = ones(nlambda,1);
xs = Xtrain*A;
xt = Xtest*A;
modelo = fitcknn(xs,Ytrain);
yhat = predict(modelo,xt);
performance(1) = mean(yhat~=Ytest);

for j=2:nlambda,
    opts.x0 = vecAt;
    [vecAt, funVal3, ValueL3]= glLeastR(X, Y, lambdas(j), opts);
    A = reshape(vecAt,d,p); A = A';
    xs = Xtrain*A;
    xt = Xtest*A;
    modelo = fitcknn(xs,Ytrain);
    yhat = predict(modelo,xt);
    performance(j) = mean(yhat~=Ytest);
end
[~,idx] = min(performance);
lambda_opt = lambdas(idx(1));
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



function out = sqrm(A)
[V,d] = eig(A);
out = V*diag(sqrt(d))*V';
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
