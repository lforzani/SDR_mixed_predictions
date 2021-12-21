% Function to compare the true parameter and the fitted parameter from
% the separate logistic regression.
function[sens spec] = ParCompare_L1idv(gamma, gamma_mat, option)
p = size(gamma,1)-1;
q = size(gamma,2);
L = size(gamma_mat,4);
sens = zeros(L,1);
spec = zeros(L,1);

if strcmp(option,'max')
    for l = 1:L
        gamma_temp = gamma_mat(:,:,:,l);
        ind = zeros(p+1,q,q);
        for j = 1:q
            for k = j:q
                    ind(:,k,j) = abs(gamma_temp(:,k,j))+abs(gamma_temp(:,j,k));
                    ind(:,j,k) = ind(:,k,j);
            end
        end
        [sens(l) spec(l)] = ParCompare_L1(gamma, ind);
    end
elseif strcmp(option, 'min')
    for l = 1:L
        gamma_temp = gamma_mat(:,:,:,l);
        ind = zeros(p+1,q,q);
        for j = 1:q
            for k = j:q
                    ind(:,k,j) = abs(gamma_temp(:,k,j)).*abs(gamma_temp(:,j,k));
                    ind(:,j,k) = ind(:,k,j);
            end
        end
        [sens(l) spec(l)] = ParCompare_L1(gamma, ind);
    end
else 
    sprintf('no specified option')
end
end
    
                    
                
   
