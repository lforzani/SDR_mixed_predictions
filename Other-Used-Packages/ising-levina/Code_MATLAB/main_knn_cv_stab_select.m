
%% knngraph Study
n = 200;
q = 10;
p = 20;

Prop = [0.2 0.5 0.8];



nfold = 5; %% CV fold number 
pthres = linspace(0,1,51);
reps = 100; %% number of repetitions in stability selection
K = 20; %% number of repetitions of experiments
for i = 4:4
        for j = 3:3
            prop = Prop(j);
        i
        j
        name = sprintf('knn_prop=%.1f_k=%d.mat',prop, i);
        load(name); %% USing the same lambda as in joint fitting
        name = sprintf('knnIdv_prop=%.1f_k=%d.mat',prop, i);
        load(name); %% USing the same lambda as in joint fitting        
        gamma_max = gamma_mat_idv;
        L = 30;  %% length of Lambda
        lambda_max = 1;
        lambda = lambda_max*0.5.^linspace(0,15,L);
        %data_y = data_y(:,:,1:3);
        %data_x = data_x(:,:,1:3);
        lambda_ind_jnt_cv = cv_select(data_y, data_x, lambda, nfold, 'joint');
        lambda_ind_sep_cv = cv_select(data_y, data_x, lambda, nfold, 'separate');
        select_prob_jnt = stability_sim_select(data_y, data_x, lambda,'joint',reps);
        select_prob_sep = stability_sim_select(data_y, data_x, lambda,'separate',reps);
       
        %% Obtain the average Sens and Spec correponding to the optimal lambda
        % CV
        SensBest_jnt_cv = mean(sens(sub2ind(size(sens), lambda_ind_jnt_cv', 1:K)));
        SpecBest_jnt_cv = mean(spec(sub2ind(size(sens), lambda_ind_jnt_cv', 1:K)));
        SensBest_sep_cv = mean(sens_max(sub2ind(size(sens_max),lambda_ind_sep_cv', 1:K)));
        SpecBest_sep_cv = mean(spec_max(sub2ind(size(sens_max), lambda_ind_sep_cv', 1:K)));
        
        % StabilitySelect
        T = length(pthres);
        sens_jnt_ss = zeros(T,1);
        spec_jnt_ss = zeros(T,1);
        sens_sep_ss = zeros(T,1);
        spec_sep_ss = zeros(T,1);

        for l = 1:length(pthres)
            gamma_ss_jnt = (select_prob_jnt>=pthres(l));
            gamma_ss_sep = (select_prob_sep>=pthres(l));
            gamma_ss_sep = select_gamma(gamma_ss_sep);
            [tmp tmp1] = ParCompare_L1(gamma, gamma_ss_jnt);
            sens_jnt_ss(l) = mean(tmp); 
            spec_jnt_ss(l) = mean(tmp1);
            [tmp tmp1] = ParCompare_L1(gamma, gamma_ss_sep);
            sens_sep_ss(l) = mean(tmp);
            spec_sep_ss(l) = mean(tmp1);
        end

        %% Save the files
        name = sprintf('knnCVSelect_prop=%.1f_k=%d.mat',prop, i);        
        save(name,  'SensBest_jnt_cv','SpecBest_jnt_cv','SensBest_sep_cv','SpecBest_sep_cv',...
                    'sens_jnt_ss','spec_jnt_ss','sens_sep_ss','spec_sep_ss',...
                    'lambda_ind_jnt_cv', 'lambda_ind_sep_cv', 'select_prob_jnt',...
                    'select_prob_sep');                                 
        end
end

        
        
        





