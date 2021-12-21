
function[sens, spec, tp, tn] = ParCompare_L1(gamma, gamma_mat)
    p = size(gamma,1)-1;
    q = size(gamma,2);
    L = size(gamma_mat,4);
    spec = zeros(L,1);
    sens = zeros(L,1);
    gamma_temp = gamma;
    theta_temp = zeros(p,q);
   
    %% remove the non-penalized value
    for j = 1:q
        gamma_temp(:,j,j) = 0;
        theta_temp(:,j) = gamma(2:p+1,j,j);
    end
    total_par = (p+1)*q*(q+1)/2-q;
    total_pos = sum(sum(sum(gamma_temp~=0)))/2+sum(sum(theta_temp~=0));
    total_neg = total_par-total_pos;
    
    
    for l = 1:L
        tpos = 0;
        tneg = 0;
    for j = 1:q
        for k = j:q
            if (j==k)
                for r = 2:p+1
                    if(gamma(r,j,k)~=0 && gamma_mat(r,j,k,l)~=0)
                        tpos = tpos+1;
                    elseif(gamma(r,j,k)==0 && gamma_mat(r,j,k,l)==0)
                        tneg = tneg+1;
                    end
                end
            else
                for r = 1:p+1
                    if(gamma(r,j,k)~=0 && gamma_mat(r,j,k,l)~=0)
                        tpos = tpos+1;
                    elseif(gamma(r,j,k)==0 && gamma_mat(r,j,k,l)==0)
                        tneg = tneg+1;
                    end
                end
            end
        end
    end
    sens(l) = tpos/total_pos;
    spec(l) = tneg/total_neg;
    end
    
               
end