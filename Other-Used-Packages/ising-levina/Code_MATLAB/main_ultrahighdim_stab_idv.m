
%% ISING Model with Covariates + L1 penalty 
n = 100;
q = 200;
p = 100;
val = 8;
L = 20;
lambda_max = 1;
lambda = lambda_max*0.5.^linspace(0,7,L);
K = 10;   %% replications
reps = 50;
pthres = linspace(0,1,51);

% name = sprintf('data_highdim_powerlaw_q%d_p%d_n%d_val%d.mat', q, p, n, val);
% load(name);
% 
% n_sub = n/2;
% %select_prob_jnt = zeros(p+1,q,q,L,K); %% record the selection freq of parameters.
% matlabpool open;
% for k = 1:4
% parfor r = 1:reps
%     sprintf('stability select reps%d round%d', k, r)
%     rng('default');
%     rng(r);
%     ind = randsample(n,n_sub);
%     ytrain = data_y(ind,:,k);
%     xtrain = data_x(ind,:,k);
%     gamma_temp = SolPath(ytrain, xtrain, lambda, 'separate');
%     name = sprintf('data_ultrahighdimIdv_q%d_p%d_n%d_val%d_reps%d_round%d.mat', q, p, n, val,k,r)
%     iSaveGamma(name, gamma_temp)
%     iClear(gamma_temp);
%     %select_prob_jnt(:,:,:,:,k) = select_prob_jnt(:,:,:,:,k)+(gamma_temp~=0);
% end
% end
% matlabpool close;


%% Organizing the results
select_prob_sepmax = zeros(p+1,q,q, K); %% record the selection freq of parameters.
select_prob_sepmin = zeros(p+1,q,q, K); %% record the selection freq of parameters.

for k = 1:K
    k
    select_max = zeros(p+1,q,q,L);
    select_min = zeros(p+1,q,q,L);
for r = 1:reps
    r
    name = sprintf('Simulations/data_ultrahighdimIdv_q%d_p%d_n%d_val%d_reps%d_round%d.mat', q, p, n, val,k,r);
    load(name);
    %size(gamma)
    [gamma_max gamma_min] = select_gamma(gamma);
    select_max = select_max+(gamma_max~=0);
    select_min = select_min+(gamma_min~=0);
end
select_max = select_max/reps;
select_max = max(select_max, [], 4);
select_prob_sepmax(:,:,:,k) = select_max;
select_min = select_min/reps;
select_min = max(select_min, [], 4);
select_prob_sepmin(:,:,:,k) = select_min;
end

clearvars -except select_prob_sepmin select_prob_sepmax;
pthres = linspace(0,1,51);
n = 100;
q = 200;
p = 100;
val = 8;
K = 10;
name = sprintf('Simulations/data_highdim_powerlaw_q%d_p%d_n%d_val%d.mat', q, p, n, val);
load(name, 'gamma');


		% StabilitySelect
        T = length(pthres);
        sens_sepmax_ss = zeros(T,K);
        spec_sepmax_ss = zeros(T,K);
        sens_sepmin_ss = zeros(T,K);
        spec_sepmin_ss = zeros(T,K);

        for l = 1:length(pthres)
            l
            gamma_tmp = (select_prob_sepmax>pthres(l));
            [sens_sepmax_ss(l,:) spec_sepmax_ss(l,:)] = ParCompare_L1(gamma, gamma_tmp);
            gamma_tmp = (select_prob_sepmin>pthres(l));
            [sens_sepmin_ss(l,:) spec_sepmin_ss(l,:)] = ParCompare_L1(gamma, gamma_tmp);
        end


        %% Save the files
		name = sprintf('Simulations/Results_UltraHighdimIdvStab_q%d_p%d_n%d_val%d.mat', q, p, n, val);
        save(name,  'select_prob_sepmax','select_prob_sepmin', 'sens_sepmax_ss', 'spec_sepmax_ss', 'sens_sepmin_ss', 'spec_sepmin_ss');   