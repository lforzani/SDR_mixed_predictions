
%% Sparsity Study
n = 200;
q = 10;
p = 20;

val = 4;
Prop = [0.2 0.5 0.8];
nEdge = [10 20 30];


nfold = 5; %% CV fold number 
pthres = linspace(0,1,51);
reps = 100; %% number of repetitions in stability selection
K = 20; %% number of repetitions of experiments
for i = 1:3
        edge = nEdge(i);         
        for j = 1:3
            prop = Prop(j);
        i
        j
         name = sprintf('Simulations/Sparse_prop=%.1f_nedge=%d.mat',prop, edge);
         load(name); %% USing the same lambda as in joint fitting
%         name = sprintf('Simulations/SparseIdv_prop=%.1f_nedge=%d.mat',prop, edge);
%         load(name); %% USing the same lambda as in joint fitting        
%         gamma_max = gamma_mat_idv;
%         L = 25;  %% length of Lambda
%         lambda_max = 1;
%         lambda = lambda_max*0.5.^linspace(0,12,L);  
%         %data_y = data_y(:,:,1:3);
%         %data_x = data_x(:,:,1:3);
%         lambda_ind_jnt_cv = cv_select(data_y, data_x, lambda, nfold, 'joint');
%         lambda_ind_sep_cv = cv_select(data_y, data_x, lambda, nfold, 'separate');
%         select_prob_jnt = stability_sim_select(data_y, data_x, lambda,'joint',reps);
%         select_prob_sep = stability_sim_select(data_y, data_x, lambda,'separate',reps);
%        
%         %% Obtain the average Sens and Spec correponding to the optimal lambda
%         % CV
%         SensBest_jnt_cv = mean(sens(sub2ind(size(sens), lambda_ind_jnt_cv', 1:K)));
%         SpecBest_jnt_cv = mean(spec(sub2ind(size(sens), lambda_ind_jnt_cv', 1:K)));
%         SensBest_sep_cv = mean(sens_max(sub2ind(size(sens_max),lambda_ind_sep_cv', 1:K)));
%         SpecBest_sep_cv = mean(spec_max(sub2ind(size(sens_max), lambda_ind_sep_cv', 1:K)));
%         
%         % StabilitySelect
%         T = length(pthres);
%         sens_jnt_ss = zeros(T,1);
%         spec_jnt_ss = zeros(T,1);
%         sens_sep_ss = zeros(T,1);
%         spec_sep_ss = zeros(T,1);
% 
%         for l = 1:length(pthres)
%             gamma_ss_jnt = (select_prob_jnt>=pthres(l));
%             gamma_ss_sep = (select_prob_sep>=pthres(l));
%             gamma_ss_sep = select_gamma(gamma_ss_sep);
%             [tmp tmp1] = ParCompare_L1(gamma, gamma_ss_jnt);
%             sens_jnt_ss(l) = mean(tmp); 
%             spec_jnt_ss(l) = mean(tmp1);
%             [tmp tmp1] = ParCompare_L1(gamma, gamma_ss_sep);
%             sens_sep_ss(l) = mean(tmp);
%             spec_sep_ss(l) = mean(tmp1);
%         end

        %% Save the files
        name = sprintf('Simulations/SparseCVSelect_prop=%.1f_nedge=%d.mat',prop, edge);  
        load(name);
        bd = 1415/3;
        q_jnt = 0;
        q_sep = 0;
        ptrue = 0;
        for j = 1:q
        ptrue =   ptrue+ sum(gamma(2:p+1,j,j)~=0)+sum(sum(gamma(:,j+1:end,j)~=0));  
        q_jnt = q_jnt + sum(sum(select_prob_jnt(2:p+1, j, j,:)))/K+ sum(sum(sum(select_prob_jnt(:,(j+1):end,j,:))))/K;
        q_sep = q_sep + sum(sum(select_prob_sep(2:p+1, j, j,:)))/K+ sum(sum(sum(select_prob_sep(:,(j+1):end,j,:))))/K;
        end
        q_jnt
        q_sep
        ptrue
        p_total = p*q+(q-1)*q/2*(p+1)
        thres_jnt = q_jnt^2/p_total/bd/2+1/2
        thres_sep = q_sep^2/p_total/bd/2+1/2
        gamma_ss_jnt = (select_prob_jnt>=thres_jnt);
        gamma_ss_sep = (select_prob_sep>=thres_sep);
        [SensBest_jnt_ss SpecBest_jnt_ss] = ParCompare_L1(gamma, gamma_ss_jnt);
        SensBest_jnt_ss = mean(SensBest_jnt_ss)
        SpecBest_jnt_ss = mean(SpecBest_jnt_ss)
        [SensBest_sep_ss SpecBest_sep_ss] = ParCompare_L1(gamma, gamma_ss_sep);
        SensBest_sep_ss = mean(SensBest_sep_ss)
        SpecBest_sep_ss = mean(SpecBest_sep_ss)

        
        save(name,  'SensBest_jnt_cv','SpecBest_jnt_cv','SensBest_sep_cv','SpecBest_sep_cv',...
                    'sens_jnt_ss','spec_jnt_ss','sens_sep_ss','spec_sep_ss',...
                    'lambda_ind_jnt_cv', 'lambda_ind_sep_cv', 'select_prob_jnt',...
                    'select_prob_sep','q_jnt', 'q_sep', 'p_total','SensBest_jnt_ss','SpecBest_jnt_ss',...
                    'SensBest_sep_ss', 'SpecBest_sep_ss');                                 
        end
end

        
        
        





