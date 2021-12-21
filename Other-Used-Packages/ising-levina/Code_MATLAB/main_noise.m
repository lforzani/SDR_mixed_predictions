
%% ISING Model with Covariates + L1 penalty 
N = [200 500];
TP = [10 20 10 20];
Q = [10 10 20 20];
P = [50 200];


K = 20;
lambda_max = 1;
L = 40;
lambda = lambda_max*0.5.^linspace(0,5,L);


for i = 1:1
    for j = 1:1
    n = N(i)
    tp = TP(j)
    q = Q(j)
        for l = 1:3
            if (l==1)
                p =tp
            else
                p = P(l-1)
            end
            name = sprintf('NoiseData_n=%d_q=%d_p=%d.mat', n, q, tp);
            load(name);
            spec = zeros(L,K);
            sens = zeros(L,K);
            gamma_mat = zeros(p+1,q,q,L,K);
            dist_1 = zeros(L,K);
            dist_0 = zeros(L,K);
            for k = 1:K
                k
                %% load the data
                y = data_y(:,:,k);
                x = data_x(:,:,k);
                if (p >tp)
                x1 = mvnrnd(zeros(p-tp,1), eye(p-tp), n);
                x = [x x1];
                end
                gamma1 = [gamma; zeros(p-tp,q,q)];
                % Obtain the solution path for ising model
                gamma_mat(:,:,:,:,k) = SolPath(y, x, lambda, 'joint');
                % Obtain Sens and Spec
                [sens(:,k), spec(:,k)] = ParCompare_L1(gamma1, gamma_mat(:,:,:,:,k));
                % Obtain the L2 distance between gamma and gamma_hat over the path
                dist_1(:,k) = DistL2(gamma1, gamma_mat(:,:,:,:,k), 1);
                dist_0(:,k) = DistL2(gamma1, gamma_mat(:,:,:,:,k), 0);
            end
            name = sprintf('RefitNoise_n=%d_q=%d_tp=%d_p=%d.mat', n, q, tp, p);
            save(name, 'gamma1','gamma_mat', 'data_y', 'data_x', 'sens', 'spec','dist_1', 'dist_0');
        end
    end
end
    





    
   
        
        
        





