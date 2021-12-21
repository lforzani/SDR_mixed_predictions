function [A,B] = getABpen4bin(C,d,invV,lambda,thr)
if nargin < 5,
    thr = 5e-2;
end
if nargin < 4,
    lambda = 0;
end
[p] = size(C,1);
maxiter = 100;
[A,B] = initAB(C,d);
if lambda > 0,
    for iter = 1:maxiter
        B = updateB(A,invV,C,B,0);
        A = updateB(B',invV,C',A',lambda); A = A';
    end
    A=orth(A);
    for j=1:p,
        if norm(A(j,:)) < thr,
            A(j,:)=0;
        end
    end
    
else
    [A,B] = getABrr(C,d,invV);
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

function B = updateB(A,invV,C,B0,lambda)
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
    
    z = [0.000 lambda];
    
    [G,w] = getGroups4bin_v3(p,d);
    opts.G = G;
    opts.ind = w;
    opts.x0 = B0;

    kronn = kron(eye(r),A);
    isqrtmV = sqrtm(invV);
    X = isqrtmV*kronn; 
    Y = isqrtmV*C(:);
    [aux, funVal3, ValueL3]= overlapping_LeastR(X, Y, z, opts);

end    

B = reshape(aux,[d,r]);
end

function A = updateA(B,invV,C,A0,lambda)
% esta funcion acutaliza A utilizando el metodo iterativo
% de minimizacion cuadratica presentado en la pagina 47
% de la tesis de sabri
%Ojo que en la formula de la tesis de Sabri (47) la identidad del kron para
%el comupto de A esta mal. En vez de I_r es I_p2 (p2=m en nuestro caso).

m = size(C,1);
d = size(B,1); 
r = size(C,2);

if lambda == 0,
    kronn = kron(B, eye(m));
    aux = inv(kronn*invV*kronn')*kronn*invV*C(:);
else
    opts=[];
    opts.init=A0;        % starting from a zero point
    opts.maxIter=1000;   % maximum number of iterations
    opts.nFlag=0;       % without normalization
    opts.rFlag=0;       % the input parameter 'rho' is a ratio in (0, 1)
    opts.rsL2=0.0;     % the squared two norm term
    opts.mFlag=1;       % smooth reformulation 
    opts.lFlag=1;       % adaptive line search
    opts.tFlag=5; 
    opts.x0 = A0;

    kronn = kron(B, eye(m));
    isqrtmV = sqrtm(invV);
    X = isqrtmV*kronn'; 
    Y = isqrtmV*C(:);
    [aux, funVal3, ValueL3]= LeastR(X, Y, lambda, opts);

end    

A = reshape(aux,[m,d]);
end

function out = sqrm(A)
[v,d] = eig(A);
out = V*diag(sqrt(d))*V';
end



