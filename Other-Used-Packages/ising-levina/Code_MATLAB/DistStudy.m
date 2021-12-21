
% Simulation Study for Identifying the Non-zero Parameters
N = 100;
P = [5 20];
Val = 2;
Prop = [0.2 0.8];
pic = 0;


for i = 1:length(N)
    for j = 1:length(P)
        for l = 1:length(Val)
            for m = 1:length(Prop)
%  n = 200;
%  p = 5;
%  val = 2;
%  prop = 0.6;

 n = N(i);
 p = P(j);
 val = Val(l);
 prop = Prop(m);
 
pic = pic+1
lambda_max = 1;
lambda  = lambda_max*0.5.^[3:1:15];
L = length(lambda);
K = 10;

spec = zeros(L,K);
sens = zeros(L,K);

%% Generate parameter

adj = textread('adj1_PowerLaw.txt');
q = size(adj,2);
gamma = gamma_gen(adj, p, val, prop);
dist = zeros(L,K);

%% Repeat experiments
for k = 1:K
    i
    j
    l
    m
    k
    
        %% Generate data
        [y,x] = DataGen0(n,gamma);
       
        %% Obtain the solution path for ising model
        [gamma_mat] = SolPath(y, x, lambda);
        
        dist(:,k) = DistCompare_L1(gamma, gamma_mat)
        
       
end


%% Make dist plot
      dist_vec = matrix(dist, K*L,1);
      lambda_vec = repmat(-log(lambda), 1,K);
      plot(lambda_vec, dist_vec, 'bo');
      hold on;
      dist1 = smooth(lambda_vec, dist_vec, 1/3, 'lowess');
      plot()
      

            end
        end
    end
end




