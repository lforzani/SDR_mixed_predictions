
%% ISING Model with Covariates + L1 penalty 
n = 200;
q = 10;
p = 20;

Val = 2.^linspace(-1,6,8);
prop = 0.5;


L = 50;  %% length of Lambda
lambda_max = 1;
lambda = lambda_max*0.5.^linspace(0,25,L);
K = 20;   %% replications
    
adj = textread(sprintf('adj%d_PowerLaw.txt', q));
gamma0 = gamma_gen(adj, p, 1, prop);

for j = 1:length(Val)
        
        j
        %% Generate Parameter 
        gamma = gamma0*Val(j);
        
        %% Repeat experiments
        spec = zeros(L,K);
        sens = zeros(L,K);
        data_y = zeros(n,q,K);
        data_x = zeros(n,p,K);
        gamma_mat = zeros(p+1,q,q,L,K);
        dist_1 = zeros(L,K);
        dist_0 = zeros(L,K);
        
        for k = 1:K
        k
        %%Generate data
        [y,x] = DataGen0(n,gamma);
        data_y(:,:,k) = y;
        data_x(:,:,k) = x;

        % Obtain the solution path for ising model
        gamma_mat(:,:,:,:,k) = SolPath(y, x, lambda, 'joint');
        % Obtain Sens and Spec
        [sens(:,k), spec(:,k)] = ParCompare_L1(gamma, gamma_mat(:,:,:,:,k));
        % Obtain the L2 distance between gamma and gamma_hat over the path
        dist_1(:,k) = DistL2(gamma, gamma_mat(:,:,:,:,k), 1);
        dist_0(:,k) = DistL2(gamma, gamma_mat(:,:,:,:,k), 0);                                                                                                                                                                                  
        end
        
        %% Save the files
        name = sprintf('Value=%.1f_n=%d_q=%d_p=%d.mat',Val(j),n,q,p);
        save(name, 'gamma', 'gamma_mat', 'data_y', 'data_x', 'sens', 'spec', 'dist_1', 'dist_0');
end


        
        
        
        





