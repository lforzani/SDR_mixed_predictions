
%% Sparsity Study
n = 200;
q = 10;
p = 20;

val = 4;
Prop = [0.2 0.5 0.8];
nEdge = [10 20 30];

K = 20;   %% replications

for i = 1:3
        edge = nEdge(i);         
        for j = 1:3
            prop = Prop(j);
        i
        j
        name = sprintf('Sparse_prop=%.1f_nedge=%d.mat',prop, edge);
        load(name);
        L = 50;  %% length of Lambda
        lambda_max = 1;
        lambda = lambda_max*0.5.^linspace(0,25,L);

        %% Repeat experiments
        spec_max = zeros(L,K);
        spec_min = spec_max;
        sens_max = spec_max;
        sens_min = spec_max;
        gamma_mat_idv = zeros(p+1,q,q,L,K);
        distIdv_1 = zeros(L,K);
        distIdv_0 = zeros(L,K);
        for k = 1:K
        k
        %%Generate data
        y = data_y(:,:,k);
        x = data_x(:,:,k);
        
        % Obtain the solution path for ising model
        gamma_mat_idv(:,:,:,:,k) = SolPath(y, x, lambda,'separate');
        
        % Obtain Sens and Spec
        [sens_max(:,k), spec_max(:,k)] = ParCompare_L1idv(gamma, gamma_mat_idv(:,:,:,:,k),'max');
        [sens_min(:,k), spec_min(:,k)] = ParCompare_L1idv(gamma, gamma_mat_idv(:,:,:,:,k),'min');
        
         % Obtain the L2 distance between gamma and gamma_hat over the path
        distIdv_1(:,k) = DistL2(gamma, gamma_mat_idv(:,:,:,:,k), 1);
        distIdv_0(:,k) = DistL2(gamma, gamma_mat_idv(:,:,:,:,k), 0);    
        end
        
        %% Obtain the average Sens and Spec correponding to the optimal lambda
        LambdaBest_jnt = Validtune(data_y, data_x, gamma_mat);
        LambdaBest_idv = Validtune(data_y, data_x, gamma_mat_idv);
       
        SensBest_jnt = mean(sens(sub2ind(size(sens), LambdaBest_jnt', 1:K)));
        SpecBest_jnt = mean(spec(sub2ind(size(sens), LambdaBest_jnt', 1:K)));
        SensBest_max = mean(sens_max(sub2ind(size(sens_max), LambdaBest_idv', 1:K)));
        SpecBest_max = mean(spec_max(sub2ind(size(sens_max), LambdaBest_idv', 1:K)));
       
        %% Save the files
        name = sprintf('SparseIdv_prop=%.1f_nedge=%d.mat',prop, edge);
        save(name,'gamma', 'gamma_mat_idv', 'lambda', 'sens_min', 'spec_min','sens_max', 'spec_max', ...
                    'distIdv_1', 'distIdv_0', 'SensBest_jnt','SpecBest_jnt','SensBest_max','SpecBest_max');
    end
end

        
        
        





