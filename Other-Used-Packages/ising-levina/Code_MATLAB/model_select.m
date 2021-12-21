% function to choose the tuning parameter by model selection criterion
% Arguments
% data_y: pool of y data (n,q,K)
% data_x: pool of x data (n,q,K)
% gamma_mat: pool of fitted parameters for each data and solution path...
%             (p+1, q, q, L, K)
function[lambda_opt eval] = model_select(data_y, data_x, gamma_mat, method, solution)
[n q K]= size(data_y);
L = size(gamma_mat, 4);
lambda_opt = zeros(K,1);
neglikeli = zeros(L,K);
eval = zeros(L,K);
for k = 1:K
    for l = 1:L
        neglikeli(l,k) = n* loss(data_y(:,:,k), data_x(:,:,k), gamma_mat(:,:,:,l,k), 0);
        npar = 0;
        if strcmp(solution, 'joint')
            for j = 1:q
            npar = npar + sum(gamma_mat(2:q+1,j,j,l,k)~=0);
            npar = npar + sum(sum(gamma_mat(:,[1:j-1 j+1:q],j,l,k)~=0))/2;
            end
        elseif strcmp(solution, 'separate')
            for j = 1:q
            tmp = sum(gamma_mat(2:q+1,j,j,l,k)~=0);
            tmp = tmp + sum(sum(gamma_mat(:,[1:j-1 j+1:q],j,l,k)~=0));
            npar = npar +tmp;
            if strcmp(method, 'BIC')
            eval(l,k) = eval(l,k)+ tmp*log(n);
            end
            end
            if strcmp(method, 'BIC')
                eval(l,k) = eval(l,k)+2*neglikeli(l,k);
            end
        end
        
        if strcmp(method, 'AIC')
            eval(l,k) = 2*neglikeli(l,k)+ 2*npar;
        elseif strcmp(method, 'BIC')
            if strcmp(solution, 'joint')
            eval(l,k) = 2*neglikeli(l,k)+ npar*log(n*q);
            end
        end
    end
    [~, lambda_opt(k)] = min(eval(:,k)); 
end
end

