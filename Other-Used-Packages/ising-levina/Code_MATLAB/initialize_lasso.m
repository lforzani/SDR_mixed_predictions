% Initialize the algorithm using glmnet package in MATLAB
%
% Call      : par_init = initialize_lasso(y,x,lambda)
%
% Arguments
% y         : n x q response matrix
% x         : n x p covariate matrix
% lambda    : penalty parameter
%
% Values
% par_init  : (p+1) x q x q initial parameter array

function[par_init] = initialize_lasso(y,x,lambda)

n = size(y,1);
q = size(y,2);
p = size(x,2);

for j = 1:q
    par_init(:,:,j) = zeros(p+1,q);
end

for j = 1:q
    for k = j:q
       for r = 1:p+1
           if(r==1)
               xx = ones(n,1);
           else
               xx = x(:,r-1);
           end
           
           if(k > j)
               xx = xx.*y(:,k);
           end
           
           if(r~=1 || k~=j)
           options = glmnetSet;
           options.nlambda = 1;
           options.lambda = lambda;
           
           fit = glmnet(xx,y(:,j)+1,'binomial',options);
           par_init(r,k,j) = fit.beta;
           par_init(r,j,k) = par_init(r,k,j);
           end
       end
    end
end

end

           

        
                    
                    
                    
                    
                    
                    
                    
                    
         
