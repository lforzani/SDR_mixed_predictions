
%% ISING Model with Covariates + L1 penalty 
n = 200;
q = 10;
p = 10;

Val = 2.^linspace(-1,3,5);
prop = 0.4;

a = figure;
pic = 1;    

    for j = 1:length(Val)
        
        j
        if(j == 1)
        name = sprintf('Value=%.1f_n=%d_q=%d_p=%d.mat',Val(j),n,q,p);
        else
        name = sprintf('Value=%d_n=%d_q=%d_p=%d.mat',Val(j),n,q,p);
        end
        load(name);
        K = 20;
        L = length(lambda);
        
        %% Make Roc plot
        sens_vec = reshape(sens, K*L, 1);
        spec_vec = reshape(spec, K*L, 1);

        subplot(3,5,j)
        plot(1-spec_vec,sens_vec, 'ko');
        axis([0,1,0,1]);
        name = sprintf('val = %.1f',Val(j));
        title(name);
        pic = pic+1;
    end
     
    for j = 1:length(Val)
        pic
        name = sprintf('ValueIdv_par=%1.1f_n=%d_q=%d_p=%d.mat',Val(j),n,q,p);
        load(name);
        K = 20;
        L = length(lambda);
        
        %% Make Roc plot
        sens_vec_max = reshape(sens_max, K*L, 1);
        spec_vec_max = reshape(spec_max, K*L, 1);
        sens_vec_min = reshape(sens_min, K*L, 1);
        spec_vec_min = reshape(spec_min, K*L, 1);

        subplot(3,5,pic);
        plot(1-spec_vec_max, sens_vec_max, 'bo');
        axis([0,1,0,1]);
        name = sprintf('val = %1.1f',Val(j));
        title(name);
        
        subplot(3,5,pic+5);
        plot(1-spec_vec_min, sens_vec_min, 'ro');
        axis([0,1,0,1]);
        name = sprintf('val = %1.1f',Val(j));
        title(name);

        pic = pic+1;
    end
  
    %% Plotting the Bayes Err
    figure;
    n = 200;
    q = 10;
    p = 10;
    pic = 1;


for j = 1:length(Val)
        if(j == 1)
        name = sprintf('Value=%.1f_n=%d_q=%d_p=%d.mat',Val(j),n,q,p);
        else
        name = sprintf('Value=%d_n=%d_q=%d_p=%d.mat',Val(j),n,q,p);
        end
        load(name);
        K = 20;
        L = length(lambda);
        
        emp0 = zeros(q,K);
        emp1 = zeros(q,K);
        for k = 1:K
        k;
        y = data_y(:,:,k);
        x = data_x(:,:,k);
        [emp0(:,k) emp1(:,k)] =BayesErr(y,x,gamma,0); 
        end
        
        subplot(2,5,pic)
        hold on;
        boxplot(emp0');
        axis([0.6,10.4,-0.1,0.5])
        m=mean(mean(emp0));
        plot(m*ones(q,1), 'k-')
        name = sprintf('Val=%.1f BayesErr1=%1.4f',Val(j), m);
        title(name)
        hold off;
        
        subplot(2,5,pic+5)
        hold on;
        boxplot(emp1');
        axis([0.6,10.4,-0.1,0.5])
        m=mean(mean(emp1));
        plot(m*ones(q,1), 'k-')
        name = sprintf('Val=%.1f BayesErr2=%1.4f',Val(j),m);
        title(name)
        hold off;
        pic = pic+1;
end




        
        
        
        





