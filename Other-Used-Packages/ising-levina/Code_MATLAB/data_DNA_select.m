%% Load data
y = csvread('DNA_Y.csv');
x = csvread('DNA_X.csv');

% subplot(1,2,1);
% hist(mean(y,1),20);
id_10 = find(mean(y)<=0.1);
y(:,id_10) = [];
% subplot(1,2,2);
% hist(mean(y,1),20)

x = zscore(x);
x(:,4) = [];

[n q] = size(y);
p = size(x,2);

%%%%%%%%%%%%%%%%%%%%% Stability selection on edges
K = 100;
L = 41;
lambda = 0.5.^linspace(0,8,L+1);
m = 0.5;
alpha = 0; %% indicator of standardization in lasso subsample
beta = 0; %% indicator of standardization in lasso before subsample
omega = 1; %% indicator of standardize y for lasso design matrix. 

%select_prob = stability_selection(y,x,lambda, m, K, alpha, beta, omega);

ind_gamma = zeros(p+1,q,q,L); %% record the selection of each 
for k = 1:K
    k
    name = sprintf('DataAnalysis/New_stdX-4_stdY_nonstdLasso_gamma_mat%d.mat',k);
    load(name, 'gamma');
    [~, gamma_tmp] = select_gamma(gamma);
    %gamma_tmp = select_gamma(gamma);
    ind_gamma = ind_gamma+(gamma_tmp~=0);
end
select_prob = zeros(q,q,p+1,L);
for l = 1:L
    for i = 1:p+1
    select_prob(:,:,i,l) = reshape(ind_gamma(i,:,:,l),q,q);
    end
end
select_prob = select_prob/K;

    select_max = max(select_prob(:,:,:,1:L),[],4);
    for i = 1:(p+1)
        name = sprintf('DataAnalysis/NewMin_stdX-4_stdY_nonstdLasso_SelectProb_X%d.csv',i-1);
        tmp = select_max(:,:,i);
        csvwrite(name, tmp);
    end
%     
% %%%%%%%%%%%%%%%%%%%%%% Node ranking: take maximum over lambda first
% % K = 100;
% % q = 430;
% % node_degree = zeros(430, 4);
% % for i = 1:3
% %     for k = 1:K
% %         k
% %     name = sprintf('DataAnalysis/stdX-4_stdY_nonstdLasso_gamma_mat%d.mat',k);
% %     load(name);
% %     x_adj = reshape(gamma_temp(i,:,:,:), q, q, 21);
% %     clear gamma_temp;
% %     x_ind = (x_adj~=0);
% %     x_ind_max = max(x_ind, [], 3);
% %     node_degree(:,i) = node_degree(:,i)+sum(x_ind_max,2)/100;
% %     end
% % end
% % csvwrite('DataAnalysis/node_degree.csv', node_degree);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%% Node ranking: take average rank for each lambda first 
% 
% K = 100;
% q = 430;
% L = 41;
% p = 4;
% 
% node_rank = zeros(100, L, 430, 4);
% for i = 1:4
%     for k = 1:K
%         k
%     name = sprintf('DataAnalysis/stdX-4_stdY_nonstdLasso_gamma_mat%d.mat',k);
%     load(name);
%     x_adj = reshape(gamma_temp(i,:,:,:), q, q, L);
%     clear gamma_temp;
%     tmp_ind = (x_adj~=0);
%     tmp = reshape(sum(tmp_ind,2),q,L);
%     for l= 1:L
%         tmp(:,l) = tiedrank(tmp(:,l));
%     end
%     tmp = 431-tmp;
%     node_rank(k,:,:,i) = reshape(tmp',1,L,q,1);
%     end
% end
% save('DataAnalysis/node_rank.mat', 'node_rank');
% 
% %%summary statistics of node_rank
% std_rank = zeros(L, 430,4);
% nintyfivepct_rank = zeros(L, 430, 4);
% seventyfivepct_rank = zeros(L, 430, 4);
% fivepct_rank = zeros(L, 430, 4);
% twentyfivepct_rank = zeros(L, 430, 4);
% mean_rank = zeros(L, 430, 4);
% median_rank = zeros(L, 430, 4);
% mad_rank = zeros(L,430,4);
% mad1_rank= mad_rank;
% for i = 1:4
%     for j = 1:q
%         j
%         tmp = node_rank(:,:,j,i);
%         std_rank(:,j,i) = std(tmp);
%         mean_rank(:,j,i) = mean(tmp, 1);
%         mad_rank(:,j,i) = mad(tmp);
%         mad1_rank(:,j,i) = mad(tmp,1);
%         tmp1 = quantile(tmp, [.05, 0.25 .50 0.75 .95]);
%         fivepct_rank(:,j,i) = tmp1(1,:);
%         twentyfivepct_rank(:,j,i) = tmp1(2,:);
%         median_rank(:,j,i) = tmp1(3,:);
%         seventyfivepct_rank(:,j,i) = tmp1(4,:);
%         nintyfivepct_rank(:,j,i) = tmp1(5,:);
%     end
% end
% 
% %names = {'std','mean', 'median', 'fivepct','twentyfivepct', 'seventyfivepct', 'nintyfivepct'};
% names = {'mad', 'median'};
% for j = 1:2
% name = char(names(j))
% figure;
% for i = 1:4
% subplot(2,2,i);
% eval(strcat('plot', '(',name, '_rank(:,:,i))'));
% title((strcat('X',sprintf('%d', i-1),'\_',name)));
% if (j==1)
%     axis([0 42 0 200])
% else
%     axis([0 42 0 450])
% end
% end
% snapnow;
% end
% 
% median_tmp = reshape(median_rank(41,:,:),430,4);
% mad_tmp = reshape(mad_rank(41,:,:), 430, 4);
% std_tmp = reshape(std_rank(41,:,:),430, 4);
% mad1_tmp = reshape(mad1_rank(41,:,:), 430,4);
% 
% for i = 2:2
%     [tmp tmp_id] = sort(median_tmp(:,i));
%     [tmp mad1_tmp(tmp_id, i) mad_tmp(tmp_id, i) std_tmp(tmp_id, i) ]
% end
% 
% csvwrite('median_rank.csv', median_tmp);

    
    
    
