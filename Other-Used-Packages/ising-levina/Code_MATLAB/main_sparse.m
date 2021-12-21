
%% Sparsity Study
n = 200;
q = 10;
p = 20;

val = 4;
Prop = [0.2 0.5 0.8];
nEdge = [10 20 30];


L = 50;  %% length of Lambda
lambda_max = 1;
lambda = lambda_max*0.5.^linspace(0,25,L);
K = 20;   %% replications

%% Generate adj matrix for edges[ 10 20 30 ]
for i = 1:length(nEdge)
    edge = nEdge(i);
    if (edge == 10)
        adj = textread(sprintf('adj%d_PowerLaw.txt', q));
    else
        adj = zeros(q,q);    
        cum = cumsum((q-1:-1:1));
        ind = randperm(q*(q-1)/2);
        ind = ind(1:edge);
        for ite = 1:edge
            ite;
            row = find(ind(ite) <= cum, 1,  'first');
            if (row>1)
                col = ind(ite)-cum(row-1)+row;
            else
                col = ind(ite)+row;
            end  
            adj(row, col) = 1;
        end
        adj = adj+adj';
    end       
    for j = 1:length(Prop)
        i
        j
        prop = Prop(j);
        gamma = gamma_gen(adj, p, val, prop);
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
        name = sprintf('Sparse_prop=%.1f_nedge=%d.mat',prop, edge);
        save(name, 'adj','gamma','gamma_mat', 'data_y', 'data_x', 'sens', 'spec','dist_1', 'dist_0');
    end
end

        
        
        
        





