

function[gamma_mat] = SolPathIdv(y, x, ind, lambda, alpha, beta, omega)
[n q] = size(y);
p = size(x,2);
L = length(lambda);
gamma_mat = zeros(p+1,q,q,L);
% gamma_mat_max = gamma_mat;


%% Separate LR with glmnet
ytmp = y;
if (omega==1)
ytmp = zscore(y);
end
K = zeros(n, (p+1)*q);
for i = 1:n
K(i,:) = kron(ytmp(i,:),[1 x(i,:)]);
end

for j = 1:q
    if(mod(j,100)==0)
    j
    end
    %% Create the data for the jth logistic regression.
    xtmp = zeros(n, p+(p+1)*(q-1));
    xtmp(:,1:p) = x;
    tmp = K;
    tmp(:,(j-1)*(p+1)+1:j*(p+1))=[];
    xtmp(:, p+1:end) = tmp;
    if (beta==1)
        xtmp = zscore(xtmp);
    end
        
    
    %% Options for glmnet
    options = glmnetSet;
    options.lambda = lambda;
    options.standardize = alpha;
    options.alpha = 1;
    
    %% Fit the solution path for given lambda of Yj
    fit = glmnet(xtmp(ind,:), y(ind,j)+1, 'binomial', options);
    B = fit.beta;
    A = fit.a0;
    clear fit xtmp;
    if (size(B,2)<L)
        A = [A A(end)*ones(1,L-size(B,2))];
        B = [B repmat(B(:,end),1,L-size(B,2))];
    end
    for l = 1:L
    beta  = B(:,l);
    beta0 = A(l);
    gamma_mat(:,:,j,l) = get_gamma(beta, beta0, p, q, j);
    end
end

% for l = 1:L
% gamma_mat_max(:,:,:,l) = select_gamma(gamma_mat(:,:,:,l));
% end
% clear gamma_mat;
end



