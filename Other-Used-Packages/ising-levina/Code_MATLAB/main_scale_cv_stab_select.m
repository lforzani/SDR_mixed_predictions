
%% ISING Model with Covariates + L1 penalty 
n = 200;
q = 10;
p = 20;

Val = 2.^linspace(-1,6,8);


nfold = 5; %% CV fold number 
pthres = linspace(0,1,51);
reps = 50; %% number of repetitions in stability selection
K = 20; %% number of repetitions of experiments
    for j = 1:4
        
%         j
%         name = sprintf('Value=%.1f_n=%d_q=%d_p=%d.mat',Val(j),n,q,p);
%         load(name); %% USing the same lambda as in joint fitting
%         name = sprintf('ValueIdv=%.1f_n=%d_q=%d_p=%d.mat',Val(j),n,q,p);
%         load(name);
%         L = 50;  %% length of Lambda
%         lambda_max = 1;
%         lambda = lambda_max*0.5.^linspace(0,25,L);  
%         %lambda_ind_jnt_cv = cv_select(data_y, data_x, lambda, nfold, 'joint');
%         %lambda_ind_sep_cv = cv_select(data_y, data_x, lambda, nfold, 'separate');
%         select_prob_jnt = stability_sim_select(data_y, data_x, lambda,'joint',reps);
%         select_prob_sep = stability_sim_select(data_y, data_x, lambda,'separate',reps);
%        
%         %% Obtain the average Sens and Spec correponding to the optimal lambda
% %         % CV
% %         SensBest_jnt_cv = mean(sens(sub2ind(size(sens), lambda_ind_jnt_cv', 1:K)));
% %         SpecBest_jnt_cv = mean(spec(sub2ind(size(sens), lambda_ind_jnt_cv', 1:K)));
% %         SensBest_sep_cv = mean(sens_max(sub2ind(size(sens_max),lambda_ind_sep_cv', 1:K)));
% %         SpecBest_sep_cv = mean(spec_max(sub2ind(size(sens_max), lambda_ind_sep_cv', 1:K)));
%         
%         % StabilitySelect
%         T = length(pthres);
%         sens_jnt_ss = zeros(T,1);
%         spec_jnt_ss = zeros(T,1);
%         sens_sep_ss = zeros(T,1);
%         spec_sep_ss = zeros(T,1);
% 
%         for i = 1:length(pthres)
%             gamma_ss_jnt = (select_prob_jnt>=pthres(i));
%             gamma_ss_sep = (select_prob_sep>=pthres(i));
%             gamma_ss_sep = select_gamma(gamma_ss_sep);
%             [tmp tmp1] = ParCompare_L1(gamma, gamma_ss_jnt);
%             sens_jnt_ss(i) = mean(tmp); 
%             spec_jnt_ss(i) = mean(tmp1);
%             [tmp tmp1] = ParCompare_L1(gamma, gamma_ss_sep);
%             sens_sep_ss(i) = mean(tmp);
%             spec_sep_ss(i) = mean(tmp1);
%         end

       %% Save the files
        name = sprintf('Simulations/ValueCVSelect=%.1f_n=%d_q=%d_p=%d.mat',Val(j),n,q,p);
        load(name);
        bd = 1145/3;
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

        
        
        





