function performance = getLambda4optimal_v2(C,d,p,q,datos,lambdas)
if nargin < 6,
    error('error: not enough input parameters');
end
if d > p || d > q,
    error('requested dimension is not compatible with the data');
end

nlambda = length(lambdas);

Xtrain = datos.Xtrain;
Wtrain = Xtrain(:,1:p);
Htrain = Xtrain(:,p+1:end);
Ytrain = datos.Ytrain;

Xtest = datos.Xtest;
Wtest = Xtest(:,1:p);
Htest = Xtest(:,p+1:end);
Ytest = datos.Ytest;

pp = d*(p+q*(q+1)/2);
[A,B] = initAB(C,d);

opts = initOpts(B',C',A',pp,p,q);
[X,Y] = getXY(B',C');

lambda = lambdas(1);
z = [0.6*lambda 0.4*lambda];
    
[vecAt]= overlapping_LeastR(X, Y, z, opts);

ppp = (p+q*(q+1)/2);
A = reshape(vecAt,d,ppp); A = A';
performance = ones(nlambda,1);
Acont = A(1:p,:);
Abin = A(p+1:end,:);
ws = (Wtrain)*Acont;
wt = (Wtest)*Acont;
hs = getT(Htrain)*Abin;
ht = getT(Htest)*Abin;
xs = [ws,hs];
xt = [wt,ht];

modelo = fitcknn(xs,Ytrain);
yhat = predict(modelo,xt);
performance(1) = mean(yhat~=Ytest);
    
for j=2:nlambda,
    opts.x0 = vecAt;
    lambda = lambdas(j);
    z = [0.6*lambda 0.4*lambda];

    [vecAt]= overlapping_LeastR(X, Y, z, opts);
    A = reshape(vecAt,d,ppp); A = A';
    Acont = A(1:p,:);
    Abin = A(p+1:end,:);
    ws = (Wtrain)*Acont;
    wt = (Wtest)*Acont;
    hs = getT(Htrain)*Abin;
    ht = getT(Htest)*Abin;
    xs = [ws,hs];
    xt = [wt,ht];
    modelo = fitcknn(xs,Ytrain);
    yhat = predict(modelo,xt);
    performance(j) = mean(yhat~=Ytest);
end
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

%     kronn = kron(eye(r),A);
%     X = kronn; 
%     Y = C(:);
%     [aux, funVal3, ValueL3]= overlapping_LeastR(X, Y, z, opts);

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

