
%% ISING Model with Covariates + L1 penalty 
N = [20 50];
TP = [10 20 10 20];
Q = [10 10 20 20];

val = 4;
prop = 0.5;
K = 20;

for i = 1:2
    for j = 1:4
    n = N(i)
    p = TP(j)
    q = Q(j)
    adj = textread(sprintf('adj%d_PowerLaw.txt', q));
    gamma = gamma_gen(adj, p, val, prop);
    data_y = zeros(n,q,K);
    data_x = zeros(n,p,K);
        for k = 1:K
            k
        [y,x] = DataGen0(n,gamma);
        data_y(:,:,k) = y;
        data_x(:,:,k) = x;
        end
    name = sprintf('NoiseData_n=%d_q=%d_p=%d.mat', n, q, p);
    save(name, 'gamma', 'data_x', 'data_y');   
    end
end





    
   
        
        
        





