
%% ISING Model with Covariates + L1 penalty 
n = 500;
Q = [10 10 10 20];
P = [10 20 40 10];


val = 2;
prop = 0.4;


L = 50;  %% length of Lambda
lambda_max = 1;
lambda = lambda_max*0.5.^linspace(0,25,L);
K = 20;   %% replications
a = figure;
    
    for j = 1:length(Q)
        q = Q(j)
        p = P(j)
        name = sprintf('n=%d_q=%d_p=%d.mat',n,q,p);
        load(name);
        %save(name, 'gamma', 'lambda', 'data_y', 'data_x', 'sens', 'spec');
        
        %% Generate Parameter 
%         adj = textread(sprintf('adj%d_PowerLaw.txt', q));
%         gamma = gamma_gen(adj, p, val, prop);
        
        %% Repeat experiments
        spec = zeros(L,K);
        sens = zeros(L,K);
%        data_y = zeros(n,q,K);
%        data_x = zeros(n,p,K);
        for k = 1:K
        k
        % Generate data
%         [y,x] = DataGen0(n,gamma);
%         data_y(:,:,k) = y;
%         data_x(:,:,k) = x;
        y = data_y(:,:,k);
        x = data_x(:,:,k);
        % Obtain the solution path for ising model
        [gamma_mat] = SolPath(y, x, lambda);
        % Obtain Sens and Spec
        [sens(:,k), spec(:,k)] = ParCompare_L1(gamma, gamma_mat);
        end
        
        %% Save the files
        name = sprintf('n=%d_q=%d_p=%d.mat',n,q,p);
        save(name, 'gamma', 'lambda', 'data_y', 'data_x', 'sens', 'spec');
        
        %% Make Roc plot
        sens_vec = reshape(sens, K*L, 1);
        spec_vec = reshape(spec, K*L, 1);

        subplot(2,2,j)
        plot(1-spec_vec,sens_vec, 'bo');
        axis([0,1,0,1]);
        name = sprintf('n=%d q=%d p=%d',n, q, p);
        title(name);

        saveas(a,sprintf('n=%d_q=%d_p=%d.fig',n,q,p))
        
    end


        
        
        
        
        





