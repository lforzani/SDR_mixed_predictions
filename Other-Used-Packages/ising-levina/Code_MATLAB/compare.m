% Function to compare the true parameter set with the estimated parameter
% set and return a contingency table with the counts of zero vs non-zeros
%
% Call   : [tab npar] = compare(par, gamma)
%
% Arguments
% par    : (p+1) x q x q estimated parameter array
% gamma  : (p+1) x q x q true parameter array
%
% Values
% tab    : 2 x 2 contingency table of zero vs non-zero  
% npar   : total number of parameters in the model

function[tab npar] = compare(par, gamma)

r = size(par,3);
q = size(par,2);
p = size(par,1)-1;
npar = q*(q+1)/2*p; % total number of parameters

tab = zeros(2,2);

for j = 1:r
    for k = j:q
        if (j~=k)
        for i = 1:(p+1)
            
            if (par(i,k,j) ~= 0 && gamma(i,k,j)~= 0)
                tab(1,1) = tab(1,1)+1;
            elseif (par(i,k,j) ~= 0 && gamma(i,k,j)== 0)
                tab(1,2) = tab(1,2)+1;
            elseif (par(i,k,j) == 0 && gamma(i,k,j)~= 0)
                tab(2,1) = tab(2,1)+1;
            else
                tab(2,2) = tab(2,2)+1;
            end
        end
        
        else
        for (i = 2:(p+1))
            if (par(i,k,j) ~= 0 && gamma(i,k,j)~= 0)
                tab(1,1) = tab(1,1)+1;
            elseif (par(i,k,j) ~= 0 && gamma(i,k,j)== 0)
                tab(1,2) = tab(1,2)+1;
            elseif (par(i,k,j) == 0 && gamma(i,k,j)~= 0)
                tab(2,1) = tab(2,1)+1;
            else
                tab(2,2) = tab(2,2)+1;
            end
        end
        end
            %tab
            %pause;
        end
    end
end


