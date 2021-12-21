
%% ISING Model with Covariates + L1 penalty 
n = 100;
q = 100;
p = 5;

val = 4;

L = 45;
lambda_max = 1;
lambda = lambda_max*0.5.^linspace(0,15,L);
K = 20;   %% replications
    

name = sprintf('Simulations/data_highdim_powerlaw_q%d_p%d_n%d_val%d.mat', q, p, n, val);
load(name);
        

        %% Repeat experiments
        spec_max = zeros(L,K);
        spec_min = spec_max;
        sens_max = spec_max;
        sens_min = spec_max;
        gamma_mat_idv = zeros(p+1,q,q,L,K);
        distIdv_1 = zeros(L,K);
        distIdv_0 = zeros(L,K);
        matlabpool open;
        parfor k = 1:K
        k
        %%load data
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
       matlabpool close;
        %[sens_max spec_max sens_min spec_min]
        %% Obtain the average Sens and Spec correponding to the optimal lambda
%         [LambdaBest_jnt likeli_jnt]= Validtune(data_y, data_x, gamma_mat);
%         [LambdaBest_idv likeli_idv]= Validtune(data_y, data_x, gamma_mat_idv);
%        
%         SensBest_jnt = mean(sens(sub2ind(size(sens), LambdaBest_jnt', 1:K)));
%         SpecBest_jnt = mean(spec(sub2ind(size(sens), LambdaBest_jnt', 1:K)));
%         SensBest_max = mean(sens_max(sub2ind(size(sens_max), LambdaBest_idv', 1:K)));
%         SpecBest_max = mean(spec_max(sub2ind(size(sens_max), LambdaBest_idv', 1:K)));
       
        %% Save the files
        name = sprintf('Simulations/Results_highdimIdv_q%d_p%d_n%d_val%d.mat',q,p,n,val);
        save(name,'gamma', 'gamma_mat_idv', 'lambda', 'sens_min', 'spec_min','sens_max', 'spec_max')
         %           'likeli_jnt', 'likeli_idv','distIdv_1', 'distIdv_0', 'SensBest_jnt','SpecBest_jnt','SensBest_max','SpecBest_max');
           
        
        
        





