
%% ISING Model with Covariates + L1 penalty 
n = 100;
q = 100;
p = 5;

val = 4;

L = 30;
lambda_max = 1;
lambda = lambda_max*0.5.^linspace(0,10,L);
K = 20;   %% replications
reps = 50;
pthres = linspace(0,1,51);
% 
% name = sprintf('Simulations/data_highdim_powerlaw_q%d_p%d_n%d_val%d.mat', q, p, n, val);
% load(name);
% 
% n_sub = n/2;
% %select_prob_jnt = zeros(p+1,q,q,L,K); %% record the selection freq of parameters.
% matlabpool open;
% for k = 1:10
% parfor r = 1:reps
%     sprintf('stability select reps%d round%d', k, r)
%     rng('default');
%     rng(r);
%     ind = randsample(n,n_sub);
%     ytrain = data_y(ind,:,k);
%     xtrain = data_x(ind,:,k);
%     gamma_temp = SolPath(ytrain, xtrain, lambda, 'joint');
%     name = sprintf('Simulations/data_highdim_q%d_p%d_n%d_val%d_reps%d_round%d.mat', q, p, n, val,k,r)
%     iSaveGamma(name, gamma_temp)
%     iClear(gamma_temp);
%     %select_prob_jnt(:,:,:,:,k) = select_prob_jnt(:,:,:,:,k)+(gamma_temp~=0);
% end
% end
% matlabpool close;

%% Organizing the results
select_prob_jnt = zeros(p+1,q,q,L,K); %% record the selection freq of parameters.

for k = 1:10
    k
for r = 1:reps
    name = sprintf('Simulations/data_highdim_q%d_p%d_n%d_val%d_reps%d_round%d.mat', q, p, n, val,k,r);
    load(name);
    select_prob_jnt(:,:,:,:,k) = select_prob_jnt(:,:,:,:,k)+(gamma~=0);
end
end

select_prob_jnt = select_prob_jnt/reps;
select_prob_jnt = max(select_prob_jnt, [], 4);
select_prob_jnt = reshape(select_prob_jnt, p+1, q,q,K);
		
name = sprintf('Simulations/data_highdim_powerlaw_q%d_p%d_n%d_val%d.mat', q, p, n, val);
load(name);
% StabilitySelect
        T = length(pthres);
        sens_jnt_ss = zeros(T,K);
        spec_jnt_ss = zeros(T,K);
        sens_sep_ss = zeros(T,K);
        spec_sep_ss = zeros(T,K);

        for l = 1:length(pthres)
            gamma_ss_jnt = (select_prob_jnt>pthres(l));
            [sens_jnt_ss(l,:) spec_jnt_ss(l,:)] = ParCompare_L1(gamma, gamma_ss_jnt);
        end

        %% Save the files
		name = sprintf('Simulations/Results_HighdimStab_q%d_p%d_n%d_val%d.mat', q, p, n, val);
        save(name,  'select_prob_jnt', 'sens_jnt_ss', 'spec_jnt_ss');   