% Initialization of the entire parameter array with individual L1-penalized
% logistic regression
%
% Call    : gg = initialize(y,x,lambda)
%
% Arguments
% y       : n x q response matrix
% x       : n x p covariate matrix
% lambda  : penalty parameter
%
% Values
% gg      : (p+1) x q x q initial parameter array

function[par_init] = initialize(y,x);



n = size(y,1);
q = size(y,2);
p = size(x,2);
par_init = zeros(p+1,q,q);
for k = 1:q
    odd = log(sum(y(:,k))/(n-sum(y(:,k))));
    if(0.01<=abs(odd) && abs(odd)<=100)
        par_init(1,k,k) = odd;
    elseif (abs(odd)>100)
        par_init(1,k,k) = 100;
    end
end

% 
% for j = 1:q
%     gg(:,:,j) = zeros(p+1,q);
% end
% 
% for j = 1:q
%     for k = j:q
%         if (k==j)  %% initialize beta's 
%             odd = log(sum(y(:,k))/(n-sum(y(:,k))));
%             if(0.01<=abs(odd) && abs(odd)<=1000)
%             gg(1,k,j) = odd;
%             end
%             for r = 1:p
%                 y1 = y(:,j);
%                 x1 = x(:,r);
%                 gg(r+1,k,j) = initialbeta(y1,x1,lambda);
%             end
%             
%         elseif (k>j) %% initialize gamma's
%             y1 = y(:,k);
%             y2 = y(:,j);
%             
%             b = y1'*y2/norm(y1,2)/norm(y2,2); %% update yj~yk
%             gg(1,k,j) = b;
%             gg(1,j,k) = gg(1,k,j);
%             
%             for r = 2:p+1                              %% update yj~x*yk
%                     x1 = x(:,r-1);
%                     gg(r,k,j) = initialgamma(y1,y2,x1,lambda);
%                     gg(r,j,k) = gg(r,k,j);
%             end
%             
%         end
%     end
%         
%         xp = [ones(n,1) x];     %% update intercepts
%  
%         A = xp*gg(:,:,j);
%         z = y;
%         z(:,j) = 1;
%         ita = diag(A*z');
%         
%         prob = cdprob(x,y,gg);
%         s = mean(log(prob(:,j)./(1-prob(:,j)))-ita);
%         if (abs(s)<10^3) && (abs(s)>10^-4)
%             gamma(1,j,j) = s;
%         end
% end
    
end

        
                    
                    
                    
                    
                    
                    
                    
                    
         
