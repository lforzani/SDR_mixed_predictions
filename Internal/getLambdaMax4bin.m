function lambda_max = getLambdaMax4bin(y,A,opts)

% compute AT y
if (opts.nFlag==0)
    ATy =A'*y;
elseif (opts.nFlag==1)
    ATy= A'*y - sum(y) * mu';  ATy=ATy./nu;
else
    invNu=y./nu;              ATy=A'*invNu-sum(invNu)*mu';
end

    lambda_max=max(abs(ATy));
    
end
