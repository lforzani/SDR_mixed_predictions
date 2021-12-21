
%% knngraph Study
n = 200;
q = 10;
p = 20;

Prop = [0.2 0.5 0.8];

L = 30;  %% length of Lambda
lambda_max = 1;
lambda = lambda_max*0.5.^linspace(0,15,L);
K = 20;   %% replications

for i = 2:4
    for j = 1:3
        i
        j
        prop = Prop(j);
        %% Repeat experiments
        spec_max = zeros(L,K);
        sens_max = zeros(L,K);
        spec_min = zeros(L,K);
        sens_min = zeros(L,K);
        name = sprintf('data_knngraph_k%d_prop%.1f.mat', i, prop);
        load(name);

        gamma_mat_idv = zeros(p+1,q,q,L,K);
        distIdv_1 = zeros(L,K);
        distIdv_0 = zeros(L,K);
        %%
        matlabpool open;
        parfor k = 1:K
        k
        
        % Obtain the solution path for ising model
        gamma_mat_idv(:,:,:,:,k) = SolPath(data_y(:,:,k), data_x(:,:,k), lambda,'separate');
        
        % Obtain Sens and Spec
        [sens_max(:,k), spec_max(:,k)] = ParCompare_L1idv(gamma, gamma_mat_idv(:,:,:,:,k),'max');
        [sens_min(:,k), spec_min(:,k)] = ParCompare_L1idv(gamma, gamma_mat_idv(:,:,:,:,k),'min');
        
         % Obtain the L2 distance between gamma and gamma_hat over the path
        distIdv_1(:,k) = DistL2(gamma, gamma_mat_idv(:,:,:,:,k), 1);
        distIdv_0(:,k) = DistL2(gamma, gamma_mat_idv(:,:,:,:,k), 0);    
        end
        matlabpool close;
        
        %% Obtain the average Sens and Spec correponding to the optimal lambda
        LambdaBest_jnt = Validtune(data_y, data_x, gamma_mat);
        LambdaBest_idv = Validtune(data_y, data_x, gamma_mat_idv);
       
        SensBest_jnt = mean(sens(sub2ind(size(sens), LambdaBest_jnt', 1:K)));
        SpecBest_jnt = mean(spec(sub2ind(size(sens), LambdaBest_jnt', 1:K)));
        SensBest_max = mean(sens_max(sub2ind(size(sens_max), LambdaBest_idv', 1:K)));
        SpecBest_max = mean(spec_max(sub2ind(size(sens_max), LambdaBest_idv', 1:K)));
       
          %% Save the files
        name = sprintf('knnIdv_prop=%.1f_k=%d.mat',prop, i);
        save(name,'gamma', 'gamma_mat_idv', 'lambda', 'sens_min', 'spec_min','sens_max', 'spec_max', ...
                    'distIdv_1', 'distIdv_0', 'SensBest_jnt','SpecBest_jnt','SensBest_max','SpecBest_max');    
    end
end

        
        
        





