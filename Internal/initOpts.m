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