
%% knngraph Study
n = 200;
q = 10;
p = 20;

Prop = [0.2 0.5 0.8];

L = 30;  %% length of Lambda
lambda_max = 1;
lambda = lambda_max*0.5.^linspace(0,15,L);
K = 20;   %% replications


for i = 2:4
    for j = 1:3
        i
        j
        prop = Prop(j);
        name = sprintf('knn_prop=%.1f_k=%d.mat',prop, i);
        load(name); %% USing the same lambda as in joint fitting
        name = sprintf('knnIdv_prop=%.1f_k=%d.mat',prop, i);
        load(name); %% USing the same lambda as in joint fitting        
        gamma_max = gamma_mat_idv;
        %% Repeat experiments
        for k = 1:K
        [gamma_max(:,:,:,:,k) gamma_min] = select_gamma(gamma_mat_idv(:,:,:,:,k));
        end
        [lambda_ind_jnt_valid likeli_sep]= Validtune(data_y, data_x, gamma_max);
        [lambda_ind_sep_valid likeli_jnt]= Validtune(data_y, data_x, gamma_mat);
        [lambda_ind_sep_aic eval_sep_aic] = model_select(data_y, data_x, gamma_max, 'AIC', 'separate');
        [lambda_ind_jnt_aic eval_jnt_aic] = model_select(data_y, data_x, gamma_mat, 'AIC', 'joint');
        [lambda_ind_sep_bic eval_sep_bic] = model_select(data_y, data_x, gamma_max, 'BIC', 'separate');
        [lambda_ind_jnt_bic eval_jnt_bic] = model_select(data_y, data_x, gamma_mat, 'BIC', 'joint');
        
        %% Obtain the average Sens and Spec correponding to the optimal lambda
       
        SensBest_jnt_valid = mean(sens(sub2ind(size(sens), lambda_ind_jnt_valid', 1:K)));
        SpecBest_jnt_valid = mean(spec(sub2ind(size(sens), lambda_ind_jnt_valid', 1:K)));
        SensBest_sep_valid = mean(sens_max(sub2ind(size(sens_max),lambda_ind_sep_valid', 1:K)));
        SpecBest_sep_valid = mean(spec_max(sub2ind(size(sens_max), lambda_ind_sep_valid', 1:K)));
        
        SensBest_jnt_aic = mean(sens(sub2ind(size(sens), lambda_ind_jnt_aic', 1:K)));
        SpecBest_jnt_aic = mean(spec(sub2ind(size(sens), lambda_ind_jnt_aic', 1:K)));
        SensBest_sep_aic = mean(sens_max(sub2ind(size(sens_max),lambda_ind_sep_aic', 1:K)));
        SpecBest_sep_aic = mean(spec_max(sub2ind(size(sens_max), lambda_ind_sep_aic', 1:K)));
        
        SensBest_jnt_bic = mean(sens(sub2ind(size(sens), lambda_ind_jnt_bic', 1:K)));
        SpecBest_jnt_bic = mean(spec(sub2ind(size(sens), lambda_ind_jnt_bic', 1:K)));
        SensBest_sep_bic = mean(sens_max(sub2ind(size(sens_max),lambda_ind_sep_bic', 1:K)));
        SpecBest_sep_bic = mean(spec_max(sub2ind(size(sens_max), lambda_ind_sep_bic', 1:K)));
       
        %% Save the files
        name = sprintf('knnModelSelect_prop=%.1f_k=%d.mat',prop, i);
        save(name,  'SensBest_jnt_valid','SpecBest_jnt_valid','SensBest_sep_valid','SpecBest_sep_valid',...
                    'SensBest_jnt_aic','SpecBest_jnt_aic','SensBest_sep_aic','SpecBest_sep_aic',...
                    'SensBest_jnt_bic','SpecBest_jnt_bic','SensBest_sep_bic','SpecBest_sep_bic',...
                    'lambda_ind_jnt_valid', 'lambda_ind_sep_valid', 'lambda_ind_jnt_aic',...
                    'lambda_ind_sep_aic', 'lambda_ind_jnt_bic', 'lambda_ind_jnt_aic',...
                    'likeli_jnt', 'likeli_sep', 'eval_sep_aic', 'eval_sep_bic', 'eval_jnt_aic', 'eval_jnt_bic');
                                 
        end
end


        
        
        





