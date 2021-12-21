% Function to fit the entire parameter array using the Coordinate Shooting 
% Algorithm on Separate Logistic Regression 
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

function[par_new] = CoordinateShooting_idv(y,x,lambda,par1)
%tic;

n = size(y,1);
q = size(y,2);
p = size(x,2);
x1 = [ones(n,1) x];

thres = 10^-3;
delta = 10^-3;
delta1 = 10^-3;
thres1 = 500;
max_ite = 10;

%% Initialization with given parameters
par_old = par1;
par_new = par_old;
prob_temp = cdprob(x,y,par_new);
w_temp = prob_temp.*(1-prob_temp);

for j = 1:q
    if (mod(j, 20)==0)
        j
    end
    dist = 1;
    loss_diff = 1;
    loss_new = 1;
    ite = 0;
    while(dist > delta && abs(loss_diff) > delta1 && ite <=max_ite && dist<=thres1)
        ite = ite+1;
        par_old(:,:,j) = par_new(:,:,j);
        for k = 1:q
            if(k==j)
                %% update theta_j0
                 u2 = sum(w_temp(:,j)); %% update beta(1,j)
                 u1 = y(:,j)-prob_temp(:,j);
                 temp = par_old(1,k,j)+sum(u1)/u2;
                    if (abs(u2)>thres)
                        par_new(1,k,j) = temp;
                    end
                prob_temp(:,j) = update_prob_idv(par_new(1,k,j),par_old(1,k,j),j,k,0,x,y,prob_temp(:,j));
                w_temp(:,j) = prob_temp(:,j).*(1-prob_temp(:,j));
%                  loss_new = loss_idv(y,x,j,par_new,lambda)
%                  pause;
                
                %% update theta_jp
                for r = 1:p                                           %% update beta (r,j)
                    u11 = u1'*x(:,r);
                    u2 = w_temp(:,j)'*(x(:,r).^2);
                    lass = lasso(u2*par_old(r+1,k,j)+u11,n*lambda);
                    temp = lass/u2;
                    if ( lass == 0 || abs(temp)<= thres)
                        par_new(r+1,k,j) = 0;
                    else
                        par_new(r+1,k,j) = temp;
                    end
                prob_temp(:,j) = update_prob_idv(par_new(r+1,k,j),par_old(r+1,k,j),j,k,r,x,y,prob_temp(:,j));
                w_temp(:,j) = prob_temp(:,j).*(1-prob_temp(:,j));
%                  loss_new = loss_idv(y,x,j,par_new,lambda)
%                  pause;
                end
            else
                %% update gamma_jkp
                
                u1 = y(:,k).*(y(:,j)-prob_temp(:,j));%n_dim
                u2 = w_temp(:,j).*y(:,k); %n_dim
                for r = 1:p+1
                   t1 = u1.*x1(:,r);
                   t2 = u2.*(x1(:,r).^2);
                   lass = lasso(sum(t1+t2*par_old(r,k,j)),n*lambda);
                   temp = lass/sum(t2);
                if( lass == 0 || abs(temp)<=thres )
                    par_new(r,k,j) = 0 ;
                else
                    par_new(r,k,j) = temp;
                end
                prob_temp(:,j) = update_prob_idv(par_new(r,k,j),par_old(r,k,j),j,k,r-1,x,y,prob_temp(:,j));
                w_temp(:,j) = prob_temp(:,j).*(1-prob_temp(:,j));
%                  loss_new = loss_idv(y,x,j,par_new, lambda)
%                  pause;
                end
            end
        end
            
            dist = sum(sum(abs(par_new(:,:,j)-par_old(:,:,j))))/(q*(p+1));
            loss_old = loss_new;
            loss_new = loss_idv(y,x,j,par_new,lambda);
            loss_diff = loss_new - loss_old;
            %pause;
    end
    j;
    %pause;
    if (isnan(dist)==1 || isnan(loss_diff)==1 || isinf(dist)||abs(loss_diff)>10 || dist > 10)
        sprintf('lambda too small')
%         j
        ite;
        dist;
        loss_diff;
        par_new(:,:,j) = par1(:,:,j);
    elseif(ite >max_ite)
        sprintf('maximum iterations reached')
        dist;
        loss_diff;
    end
end
%toc;
end
                    
                   
                 
                
              
            
         

               
               
             
                    
        
    
    
  
  



