
function[gamma_max gamma_min] = select_gamma(gamma)
q = size(gamma,2);
p = size(gamma,1);
gamma_max = gamma;
gamma_min = gamma;
L = size(gamma,4);
for l = 1:L
for j = 1:q-1
    for k = j+1:q
        for i = 1:p
            if(abs(gamma(i,k,j,l)>abs(gamma(i,j,k,l))))
                gamma_min(i,k,j,l) = gamma(i,j,k,l);
                gamma_max(i,j,k,l) = gamma(i,k,j,l);
            else
                gamma_max(i,k,j,l) = gamma(i,j,k,l);
                gamma_min(i,j,k,l) = gamma(i,k,j,l);
            end
        end
    end
end
end

end
        


