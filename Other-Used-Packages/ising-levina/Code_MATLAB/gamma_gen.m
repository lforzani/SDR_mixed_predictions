
function[gamma] = gamma_gen(adj, p, val, prop)
q = size(adj,1);
gamma = zeros(p+1,q,q);
for j = 1:q
    for k = j:q
    if (k == j)
        gamma(:,j,j) = val*binornd(1, prop, p+1,1).*(2*binornd(1,0.5,p+1,1)-1);
        gamma(1,j,j) = val*(2*binornd(1,0.5)-1);
%         gamma(:,j,j) = val*binornd(1, prop, p+1,1);
%         gamma(1,j,j) = val;
    elseif (adj(j,k) == 1)
        gamma(:,k,j) = val*binornd(1, prop, p+1,1).*(2*binornd(1,0.5,p+1,1)-1);
        %gamma(:,k,j) = val*binornd(1, prop, p+1,1);
        gamma(:,j,k) = gamma(:,k,j);
    end
    end
end
end

        
        
        

