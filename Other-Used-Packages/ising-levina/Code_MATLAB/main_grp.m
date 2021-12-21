%adj = textread('adj5_PowerLaw.txt');

q = 5;
adj = ones(q,q);

p = 40;
mm = [10 20 30];
val = 2;
prop = 0.2;
n = 400;

for k = 1:3
    m = mm(k);
    gamma = gammagen_grp(adj, p, val, m, prop);
    [y,x] = DataGen0(n,gamma);

par_init = initialize(y,x,lambda1(end));

lambda1_max = 1;
lambda2_max = 1;
lambda1 = lambda1_max* 0.5.^(0:1:15);
lambda2 = lambda2_max* 0.5.^(0:1:15); 
L1 = length(lambda1);
L2 = length(lambda2);
lambda = zeros(L1*L2,2);
gamma_mat = zeros(p+1, q, q, L1*L2);
cnt = 1;
for i = 1:length(lambda1)
    for j = 1:length(lambda2)
        i
        j
        lambda(cnt,:) = [lambda1(i), lambda2(j)];
        gamma_mat(:,:,:,cnt) = GroupLasso(y,x,lambda1(i), lambda2(j), par_init);
        cnt = cnt + 1;
    end
end

[sens_grp spec_grp sens_idv spec_idv] = ParCompare_grp(gamma, gamma_mat, L1, L2);

figure;
for i  = 1:L1
    subplot(4,4,i)
    plot(1-spec_grp(:,i), sens_grp(:,i), 'o-');
    title(sprintf('lambda_1 = %d', lambda1(i)));
    xlim([0 1]);
    ylim([0 1]);
end

figure;
for j = 1:L2
    subplot(4,4,j)
    plot(1-spec_idv(j,:), sens_idv(j,:), 'o-');
    title(sprintf('lambda_2 = %d', lambda2(j)));
    xlim([0 1]);
    ylim([0 1]);
end
end
        
