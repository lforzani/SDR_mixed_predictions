
%% Simulation for ultra high dimensional 
n = 100;
q = 200;
p = 100;

val = 8;

L = 25;
lambda_max = 1;
lambda = lambda_max*0.5.^linspace(0,7,L);
K = 8;   %% replications
    
        
        %% Repeat experiments
        spec = zeros(L,K);
        sens = zeros(L,K);
        data_y = zeros(n,q,K);
        data_x = zeros(n,p,K);
        gamma_mat = zeros(p+1,q,q,L,K);
        dist_1 = zeros(L,K);
        dist_0 = zeros(L,K);
        name = sprintf('Simulations/data_highdim_powerlaw_q%d_p%d_n%d_val%d.mat', q, p, n, val);
        load(name);
        matlabpool open;
        parfor k = 1:K
            k      
        y = data_y(:,:,k);
        x = data_x(:,:,k);
        % Obtain the solution path for ising model
        gamma_mat(:,:,:,:,k) = SolPath(y, x, lambda, 'joint');
        % Obtain Sens and Spec
        [sens(:,k), spec(:,k)] = ParCompare_L1(gamma, gamma_mat(:,:,:,:,k));
        % Obtain the L2 distance between gamma and gamma_hat over the path
        dist_1(:,k) = DistL2(gamma, gamma_mat(:,:,:,:,k), 1);
        dist_0(:,k) = DistL2(gamma, gamma_mat(:,:,:,:,k), 0);                                                                                                                                                                                  
        end
        matlabpool close;
  
        %% Save the files
        name = sprintf('Simulations/Results_ultra_highdim_q%d_p%d_n%d_val%d.mat',q,p,n,val);
        save(name, 'gamma', 'gamma_mat', 'data_y', 'data_x', 'sens', 'spec', 'dist_1', 'dist_0');

        
        
        
        





