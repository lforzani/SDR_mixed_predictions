% Function to fit the entire parameter array using the Coordinate Shooting 
% Algorithm
%
% Call    : par_new = CoordinateShooting(y,x,lambda1,lambda2,par_init)
%
% Arguments
% y       : n x q response matrix
% x       : n x p covariate matrix
% lambda1 : penalty parameter for L1-penalty
% lambda2 : penalty parameter for group-lasso penalty
% par_init  : (p+1) x q x q parameter array, the initial value
%
% Values
% par_new  : (p+1) x q x q fitted parameter array

function[par_new] = GroupLasso(y,x,lambda1, lambda2 ,par_init)
tic;

n = size(y,1);
q = size(y,2);
p = size(x,2);
x1 = [ones(n,1) x];

thres = 10^-3;
thres_max = 100;
npar = q*(q+1)*(p+1)/2;
delta = 10^-3;
dist = 10;

%% Initialization with given parameters
par_old = par_init;


par_new = par_old;
prob_temp = cdprob(x,y,par_new);
w_temp = prob_temp.*(1-prob_temp);
run = 0;
los  = loss_grp(y,x,par_new, lambda1, lambda2);
los_old = los-1;
df = q*(q-1)/2;
ite = 0;
% while((dist >= delta) && (abs(los_old - los)>= delta))
while(dist >= delta && ite <=1000)
    ite = ite+1;
dist = 0;
par_old = par_new;

square = zeros(p+1,1);
for m = 1:q
    square = square + sum(par_old(:,m+1:q,m).^2 ,2);
end
square = (square.^0.5).*(square > thres);
%pause;


for j=1:q
        for k= j:q    
            if (k>j)  
               u1 = y(:,k).*(y(:,j)-prob_temp(:,j))+ y(:,j).*(y(:,k)-prob_temp(:,k));%n_dim
               u2 = w_temp(:,j).*y(:,k)+w_temp(:,k).*y(:,j); %n_dim
               
               r = 1;  % udpate gamma(j,k,0)
               t1 = u1.*x1(:,r);
               t2 = u2.*(x1(:,r).^2);
               temp = lasso(sum(t1+t2*par_old(r,k,j)),n*lambda1)/sum(t2);
               if( abs(temp)<=thres )
                    par_new(r,k,j) = 0 ;
               else
                    par_new(r,k,j) = temp;
               end
               par_new(r,j,k) = par_new(r,k,j);
               prob_temp = update_prob(par_new(r,k,j),par_old(r,k,j),j,k,r-1,x,y,prob_temp);
               w_temp = update_w(w_temp,prob_temp,j,k);
           
               
               %% update gamma(j,k,:)
               for r = 2:p+1
                   if(square(r)>0)
                   t1 = -u1'*x1(:,r)*square(r)*df;
                   t2 = u2'*(x1(:,r).^2)*square(r)*df;
                   temp = par_old(r,j,k)- (t1+n*lambda2*par_old(r,j,k))/(t2+n*lambda2); %% Newton-Raphson Update
                   par_new(r,k,j) = temp;
                   else
                       par_new(r,k,j) = 0;
                   end
                   par_new(r,j,k) = par_new(r,k,j);
                   prob_temp = update_prob(par_new(r,k,j),par_old(r,k,j),j,k,r-1,x,y,prob_temp);
                   w_temp = update_w(w_temp,prob_temp,j,k);
               end
                       
             elseif (k == j)
                    
               ss = sum(w_temp(:,j)); %% update beta(1,j)
               temp = par_old(1,k,j)+sum(y(:,j)-prob_temp(:,j))/ss;
               if(abs(temp)<=thres_max)
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
                    temp = lasso(s+t*par_old(r+1,k,j),n*lambda1)/t;
                    if (abs(temp)<= thres)   
                        par_new(r+1,k,j) = 0;
                    elseif (abs(temp)>thres && abs(temp)<=thres_max)
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
       dist = dist+sum(sum(A(:,j:q,j)));
   end
   %dist
   los = loss_grp(y,x,par_new,lambda1, lambda2);
   %pause;
   %run = run+1
   
   if(mod(run,200)==0)
       los_old = los;
       los = loss_grp(y,x,par_new,lambda1, lambda2);
       if (los - los_old >=10^-2)
       par_new = par_old;
       break;
       end
       %pause;
   end

   
end
%run
toc
end