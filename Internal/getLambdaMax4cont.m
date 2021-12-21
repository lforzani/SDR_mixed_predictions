function lambda_max = getLambdaMax4cont(y,A,opts)

% compute AT y
if (opts.nFlag==0)
    ATy =A'*y;
elseif (opts.nFlag==1)
    ATy= A'*y - sum(y) * mu';  ATy=ATy./nu;
else
    invNu=y./nu;              ATy=A'*invNu-sum(invNu)*mu';
end

if (~isfield(opts,'q'))
    q=2; opts.q=2;
else
    q=opts.q;
    if (q<1)
        error('\n q should be larger than 1');
    end
end
  
    if q==1
        q_bar=Inf;
    elseif q>=1e6
        q_bar=1;
    else
        q_bar=q/(q-1);
    end
    
    % compute the norm of ATy corresponding to each group
    ind=opts.ind;
    k=length(ind)-1; % the number of groups

    norm_ATy=zeros(k,1);
    for i=1:k
        norm_ATy(i,1)=norm(  ATy( (1+ind(i)):ind(i+1) ), q_bar );
    end
    
    % incorporate the gWeight
% gWeight: the weigtht for each group
if (isfield(opts,'gWeight'))
    gWeight=opts.gWeight;
    if (size(gWeight,1)~=k)
        error('\n opts.gWeight should a %d x 1 vector',k);
    end
    
    if (min(gWeight)<=0)
        error('\n .gWeight should be positive');
    end
else
    gWeight=ones(k,1);
end
    norm_ATy=norm_ATy./gWeight;
    
    % compute lambda_max
    lambda_max=max(norm_ATy);
    
end
