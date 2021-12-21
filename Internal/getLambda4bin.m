function [lambda_opt,performance] = getLambda4bin(C,d,datos,nlambda)
if nargin < 4,
    nlambda = 20;
end

[p] = size(C,1);
[A,B] = initAB(C,d);
opts = initOpts4bin(B',C',A');
[X,Y] = getXY(B',C');

if length(nlambda)==1,
    lambda_max = getLambdaMax4bin(Y,X,opts);
    z1 = logspace(-1,-5,nlambda);
    lambdas = z1*lambda_max;
else
    lambdas = nlambda;
end

Xtrain = datos.Xtrain;
Ytrain = datos.Ytrain;
Xtest = datos.Xtest;
Ytest = datos.Ytest;
TXtrain = getT(Xtrain);
TXtest = getT(Xtest);


lambda = lambdas(1);
z = [0.6*lambda 0.4*lambda];
    
[vecAt]= overlapping_LeastR(X, Y, z, opts);
 
A = reshape(vecAt,d,p); A = A';
nlambda = length(lambdas);
performance = ones(nlambda,1);
xs = (TXtrain)*A;
xt = (TXtest)*A;
modelo = fitcknn(xs,Ytrain);
yhat = predict(modelo,xt);
performance(1) = mean(yhat~=Ytest);
    
for j=2:nlambda,
    opts.x0 = vecAt;
    lambda = lambdas(j);
    z = [0.6*lambda 0.4*lambda];

    [vecAt]= overlapping_LeastR(X, Y, z, opts);
    A = reshape(vecAt,d,p); A = A';
    xs = (TXtrain)*A;
    xt = (TXtest)*A;
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

function B = updateB(A,C,B0,lambda)
% esta funcion acutaliza B utilizando el metodo iterativo
% de minimizacion cuadratica presentado en la pagina 47
% de la tesis de sabri
[p,r] = size(C);
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
    
    z = [0.6*lambda 0.4*lambda];
    
    [G,w] = getGroups4bin(p,d);
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