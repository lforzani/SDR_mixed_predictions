
%% knngraph Study
n = 200;
q = 10;
p = 20;

Prop = [0.2 0.5 0.8];

L = 30;  %% length of Lambda
lambda_max = 1;
lambda = lambda_max*0.5.^linspace(0,15,L);
K = 20;   %% replications

%% model fitting
for i = 2:4
    for j = 1:3
        i
        j
        prop = Prop(j);
        %% Repeat experiments
        spec = zeros(L,K);
        sens = zeros(L,K);
        gamma_mat = zeros(p+1,q,q,L,K);
        dist_1 = zeros(L,K);
        dist_0 = zeros(L,K);
        name = sprintf('data_knngraph_k%d_prop%.1f.mat', i, prop);
        load(name);
        matlabpool open;
        parfor k = 1:K
        k
        % Obtain the solution path for ising model
        gamma_mat(:,:,:,:,k) = SolPath(data_y(:,:,k), data_x(:,:,k), lambda, 'joint');
        % Obtain Sens and Spec
        [sens(:,k), spec(:,k)] = ParCompare_L1(gamma, gamma_mat(:,:,:,:,k));
        % Obtain the L2 distance between gamma and gamma_hat over the path
        dist_1(:,k) = DistL2(gamma, gamma_mat(:,:,:,:,k), 1);
        dist_0(:,k) = DistL2(gamma, gamma_mat(:,:,:,:,k), 0);    
        end
        matlabpool close;
        
        %% Save the files
        name = sprintf('knn_prop=%.1f_k=%d.mat',prop, i);
        save(name, 'adj','gamma','gamma_mat', 'data_y', 'data_x', 'sens', 'spec','dist_1', 'dist_0');
    end
end

        
        
        
        





