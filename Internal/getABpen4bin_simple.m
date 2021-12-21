function [A,B] = getABpen4bin_simple(C,d,lambda,thr)
if nargin < 4,
    thr = 1e-3;
end
if nargin < 3,
    lambda = 0;
end
[p] = size(C,1);
maxiter = 100;
[A,B] = initAB(C,d);
if lambda > 0,
%     for iter = 1:maxiter
%         B = updateB(A,invV,C,B,0);
%         A = updateB(B',invV,C',A',lambda); A = A';
%     end
    A = updateB(B',C',A',lambda); A = A';
    A=orth(A);
    for j=1:p,
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



