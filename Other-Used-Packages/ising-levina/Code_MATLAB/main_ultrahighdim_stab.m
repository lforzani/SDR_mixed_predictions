% %matlab example for matlab pool
% %this example assumes the MPIEXEC parallel config
% sched = findResource('scheduler', 'type', 'mpiexec')
% set(sched, 'MpiexecFileName', '/home/software/rhel6/mpiexec/bin/mpiexec')
% set(sched, 'EnvironmentSetMethod', 'setenv')
% %use the 'sched' object when calling matlabpool
% %the size of the pool (4) should equal ppn in PBS
% matlabpool (sched, 15)
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
% 
% name = sprintf('data_highdim_powerlaw_q%d_p%d_n%d_val%d.mat', q, p, n, val);
% load(name);
% 
% n_sub = n/2;
% %select_prob_jnt = zeros(p+1,q,q,L,K); %% record the selection freq of parameters.
% %matlabpool open;
% for k = 7:10
% parfor r = 1:reps
%     sprintf('stability select reps%d round%d', k, r)
%     rng('default');
%     rng(r);
%     ind = randsample(n,n_sub);
%     ytrain = data_y(ind,:,k);
%     xtrain = data_x(ind,:,k);
%     gamma_temp = SolPath(ytrain, xtrain, lambda, 'joint');
%     name = sprintf('data_ultrahighdim_q%d_p%d_n%d_val%d_reps%d_round%d.mat', q, p, n, val,k,r)
%     iSaveGamma(name, gamma_temp)
%     iClear(gamma_temp);
%     %select_prob_jnt(:,:,:,:,k) = select_prob_jnt(:,:,:,:,k)+(gamma_temp~=0);
% end
% end
% matlabpool close;
%% Organizing the results
select_prob_jnt = zeros(p+1,q,q,L, K); %% record the selection freq of parameters.

for k = 1:10
    k
for r = 1:reps
    r
    name = sprintf('Simulations/data_ultrahighdim_q%d_p%d_n%d_val%d_reps%d_round%d.mat', q, p, n, val,k,r);
    load(name);
    select_prob_jnt(:,:,:,:,k) = select_prob_jnt(:,:,:,:,k)+(gamma~=0);
    clear gamma;
end
end

clearvars -except select_prob_jnt

name = sprintf('Simulations/data_highdim_powerlaw_q200_p100_n100_val8.mat');
load(name, 'gamma');

select_prob_jnt = max(select_prob_jnt, [], 4);
q = 200;
p = 100;
K = 10;
reps = 50;
n = 100;
val = 8;
select_prob_jnt = reshape(select_prob_jnt, p+1, q,q,K);
select_prob_jnt = select_prob_jnt/reps;

		
pthres = linspace(0,1,51);
		% StabilitySelect
        T = length(pthres);
        sens_jnt_ss = zeros(T,K);
        spec_jnt_ss = zeros(T,K);
        sens_sep_ss = zeros(T,K);
        spec_sep_ss = zeros(T,K);

        for l = 1:length(pthres)
            l
            gamma_ss_jnt = (select_prob_jnt>pthres(l));
            [sens_jnt_ss(l,:) spec_jnt_ss(l,:)] = ParCompare_L1(gamma, gamma_ss_jnt);
        end
        

        %% Save the files
		name = sprintf('Simulations/Results_UltraHighdimStab_q%d_p%d_n%d_val%d.mat', q, p, n, val);
        save(name,  'select_prob_jnt', 'sens_jnt_ss', 'spec_jnt_ss');   