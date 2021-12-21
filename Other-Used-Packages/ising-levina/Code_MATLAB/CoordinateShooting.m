% Function to fit the entire parameter array using the Coordinate Shooting 
% Algorithm
%
% Call    : par_new = CoordinateShooting(y,x,lambda,par1)
%
% Arguments
% y       : n x q response matrix
% x       : n x p covariate matrix
% lambda  : penalty parameter
% par1    : (p+1) x q x q parameter array, the initial value
%
% Values
%par_new  : (p+1) x q x q fitted parameter array

function[par_new] = CoordinateShooting(y,x,lambda,par1)
%tic;

n = size(y,1);
q = size(y,2);
p = size(x,2);
x1 = [ones(n,1) x];
npar = (p+1)*q*(q+1)/2;

thres = 10^-3;
delta = 10^-3;
delta1 = 10^-3;
dist = 10;
max_ite = 10;

%% Initialization with given parameters
par_old = par1;
par_new = par_old;
prob_temp = cdprob(x,y,par_new);
w_temp = prob_temp.*(1-prob_temp);

[los like]  = loss(y,x,par_new, lambda);
loss_diff = 1;
like_diff = 1;
ite = 0;

while((dist >= delta) && (abs(loss_diff) >=delta1) &&(ite <= max_ite))
%while((dist >= delta) && (abs(like_diff) >=delta1) &&(ite <= max_ite))
%while((dist >= delta)  &&(ite <= max_ite))
dist = 0;
par_old = par_new;
%ite
    for j=1:q
        for k= j:q    
            if (k>j)  
               u1 = y(:,k).*(y(:,j)-prob_temp(:,j))+ y(:,j).*(y(:,k)-prob_temp(:,k));%n_dim
               u2 = w_temp(:,j).*y(:,k)+w_temp(:,k).*y(:,j); %n_dim
               for r = 1:p+1
                   t1 = u1.*x1(:,r);
                   t2 = u2.*(x1(:,r).^2);
                   lass = lasso(sum(t1+t2*par_old(r,k,j)),n*lambda);
                   temp = lass/sum(t2);
                   
               if( lass == 0 || abs(temp)<=thres )
                    par_new(r,k,j) = 0 ;
               %elseif (abs(temp) > thres && abs(temp) <= thres_max)
               else
                    par_new(r,k,j) = temp;
               end
               par_new(r,j,k) = par_new(r,k,j);
               prob_temp = update_prob(par_new(r,k,j),par_old(r,k,j),j,k,r-1,x,y,prob_temp);
               w_temp = update_w(w_temp,prob_temp,j,k);
               %j
               %k
               %1
%                %los = loss(y,x,par_new,lambda) m,./
               end
              
                
           elseif (k == j)
                
               
               ss = sum(w_temp(:,j)); %% update beta(1,j)
               temp = par_old(1,k,j)+sum(y(:,j)-prob_temp(:,j))/ss;
               if (abs(ss)<=thres)
                   par_new(1,k,j) = 0;
               else
                   par_new(1,k,j) = temp;
               end
               
             
                    prob_temp = update_prob(par_new(1,k,j),par_old(1,k,j),j,k,0,x,y,prob_temp);
                    w_temp = update_w(w_temp,prob_temp,j,k);
                    %j
                    %k
                    %1
                    %los = loss(y,x,par_new,lambda)
               
               
             
               for r = 1:p                                           %% update beta (r,j)
                    s = (y(:,j)-prob_temp(:,j))'*x(:,r);
                    t = w_temp(:,j)'*(x(:,r).^2);
                    lass = lasso(s+t*par_old(r+1,k,j),n*lambda);
                    temp = lass/t;
                    if ( lass == 0 || abs(temp)<= thres)
                        par_new(r+1,k,j) = 0;
                    %elseif (abs(temp)>thres && abs(temp)<=thres_max)
                    else
                        par_new(r+1,k,j) = temp;
                    end
                    
                   
                    prob_temp = update_prob(par_new(1+r,k,j),par_old(1+r,k,j),j,k,r,x,y,prob_temp);
                    w_temp = update_w(w_temp,prob_temp,j,k);
                    %j
                    %k
                    %r+1
                    %los = loss(y,x,par_new,lambda)
               end
               %pause;
               
               
              
           end
           
       end
       
   end
   A = abs(par_new-par_old);
   
   for j = 1:q
       dist = dist+sum(sum(A(:,j:q,j)))/npar;
   end
   
    %if(mod(ite,10)==0)
       dist;
       los_old = los;
       like_old = like;
       [los like] = loss(y,x,par_new,lambda);
       loss_diff = los_old - los;
       like_diff = like_old - like;
       %pause;
    %end
    ite = ite+1;
    end

if (isnan(dist)==1 || isnan(loss_diff)==1)
    %sprintf('lambda too small')
    %par_new = par1;
end

if(ite >max_ite)
   % sprintf('maximum iterations reached')
end
%run
%toc
end
