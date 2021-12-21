%function [misclass_optimal,misclass_optimal_t,misclass_suboptimal,misclass_suboptimal_t,misclass_pca,misclass_all] = cvPFCmix(Y,X,H,type,model,nfold)

function [tabla_optimal, tabla_suboptimal, tabla_all] = cvPFCmix(Y,X,H,type,model,nfold)

if nargin<6, %ojo, cambiar al final 
    nfold=5;
end
disp('== starting CV procedure for PFCMix (this can take several minutes)===')

n = length(Y);
p = size(X,2);
q = size(H,2);
kk = q*(q+1)/2;
% 
% err1_optimal = zeros(1,nfold);
% err2_optimal = zeros(1,nfold);
% err1_suboptimal = zeros(1,nfold);
% err2_suboptimal = zeros(1,nfold);
% err1_all = zeros(1,nfold);
% err2_all = zeros(1,nfold);
% err1_nlpca = zeros(1,nfold);
% err2_nlpca = zeros(1,nfold);
  
err_optimal = zeros(1,nfold);
err_suboptimal = zeros(1,nfold);
err_optimal_t = zeros(1,nfold);
err_suboptimal_t = zeros(1,nfold);
err_all = zeros(1,nfold);
err_pca = zeros(1,nfold);

yhat_optimal = zeros(1,nfold);
yhat_suboptimal = zeros(1,nfold);
yhat_all = zeros(1,nfold);
        
blksize = floor(n/nfold);
shuf = randperm(n);
Y = Y(shuf); X = X(shuf,:); H = H(shuf,:);

% t_optimal_aux = zeros(2,2,nfold);
% t_suboptimal_aux = zeros(2,2,nfold);
% t_all_aux = zeros(2,2,nfold);

for k = 1:nfold
    k
    % set training and testing partitions
    if k < nfold
        idx = ((k-1)*blksize + 1):(k*blksize);
    else
        idx = ((k-1)*blksize + 1):n;
    end
    xts = X(idx,:); yts = Y(idx,:); hts = H(idx,:);
    xtr = X; xtr(idx,:) = []; ytr = Y; ytr(idx)=[]; htr = H; htr(idx,:)=[];
    nts = size(xts,1);
    ntr = size(xtr,1);

    % estimate dimension
    %[redu_optimal,redu_suboptimal,fycent,~,~,Tauhatraro0, Tauhatraro, Ahat, betahat, Deltahat] = EM4mixture_continua_binaria(xtr,htr,ytr,type);

    
%     [dtest1, dtest2, dtestfirst1, dtestfirst2, dtestsecond1, dtestsecond2] = dselection(Tauhatraro0,Deltahat,Ahat,Tauhatraro,betahat,fycent,redu_optimal,redu_suboptimal);
%     dim = min(dtest1, dtest2);
%     dim1 = min(dtestfirst1,dtestsecond1);
%     dim2 = min(dtestfirst2,dtestsecond2)

    dim = 1;
    dim1 = 1;
    dim2 = 1;
    % estimate the reduction
    red_est_optimal = rr4mix_optimal(xtr,htr,ytr,dim,type,'auto');
    [red_est_suboptimal1,red_est_suboptimal2] = rr4mix_suboptimal(xtr,htr,ytr,dim1,dim2,type,'auto');
    
    red_est_optimal_t = rr4mix_optimal(xtr,htr,ytr,dim,type,0);
    [red_est_suboptimal1_t,red_est_suboptimal2_t] = rr4mix_suboptimal(xtr,htr,ytr,dim1,dim2,type,0,0);
   
    if(size(red_est_suboptimal2,2) == 0)
        red_est_suboptimal2 = zeros(kk, 1);
    end
   % red_est_suboptimal = [red_est_suboptimal1; red_est_suboptimal2];

    Hraro_ts = zeros(nts,kk);
    for j = 1:nts
        auxts = hts(j,:)'*hts(j,:);
        Hraro_ts(j,:) = [diag(auxts);auxts(find(tril(ones(q,q),-1)))];
    end
    T_optimal_ts = [xts,Hraro_ts];
    T_suboptimal_ts1 = [xts,hts];
    T_suboptimal_ts2 = Hraro_ts;

    Hraro_tr = zeros(nts,kk);
    for j = 1:ntr
        auxtr = htr(j,:)'*htr(j,:);
        Hraro_tr(j,:) = [diag(auxtr);auxtr(find(tril(ones(q,q),-1)))];
    end
    T_optimal_tr = [xtr,Hraro_tr];
    T_suboptimal_tr1 = [xtr,htr];
    T_suboptimal_tr2 = [Hraro_tr];
         
    vtr = [xtr htr];
    vts = [xts hts];
    
    % proyecciones
    Proy_est_optimal_tr = T_optimal_tr*red_est_optimal;
    Proy_est_optimal_ts = T_optimal_ts*red_est_optimal;
    Proy_est_suboptimal_tr1 = T_suboptimal_tr1*red_est_suboptimal1; 
    Proy_est_suboptimal_tr2 = T_suboptimal_tr2*red_est_suboptimal2;
    Proy_est_suboptimal_ts1 = T_suboptimal_ts1*red_est_suboptimal1; 
    Proy_est_suboptimal_ts2 = T_suboptimal_ts2*red_est_suboptimal2; 
    
    Proy_est_suboptimal_tr = [Proy_est_suboptimal_tr1,Proy_est_suboptimal_tr2];
    Proy_est_suboptimal_ts = [Proy_est_suboptimal_ts1,Proy_est_suboptimal_ts2];
    
    Proy_est_optimal_t_tr = T_optimal_tr*redu_optimal;
    Proy_est_optimal_t_ts = T_optimal_ts*redu_optimal;
    Proy_est_suboptimal_t_tr1 = T_suboptimal_tr1*red_est_suboptimal1_t; 
    Proy_est_suboptimal_t_tr2 = T_suboptimal_tr2*red_est_suboptimal2_t;
    Proy_est_suboptimal_t_ts1 = T_suboptimal_ts1*red_est_suboptimal1_t; 
    Proy_est_suboptimal_t_ts2 = T_suboptimal_ts2*red_est_suboptimal2_t; 
   
    Proy_est_suboptimal_t_tr = [Proy_est_suboptimal_t_tr1,Proy_est_suboptimal_t_tr2];
    Proy_est_suboptimal_t_ts = [Proy_est_suboptimal_t_ts1,Proy_est_suboptimal_t_ts2];
    
    
    % PCA
%    [~,score] = pca(vtr);
%    pc = score(1,:);
%    Proy_pca_tr = vtr*pc';
%    Proy_pca_ts = vts*pc';
%      
    if strcmpi(model,'linear'),   % fit a linear model
        m_optimal = fitlm(Proy_est_optimal_tr,ytr);
        m_suboptimal = fitlm(Proy_est_suboptimal_tr,ytr);
        m_optimal_t = fitlm(Proy_est_optimal_t_tr,ytr);
        m_suboptimal_t = fitlm(Proy_est_suboptimal_t_tr,ytr) ;
        mm = fitlm(vtr, ytr);
        %m_pca = fitlm(Proy_pca_tr, ytr);
        yhat_optimal = predict(m_optimal,Proy_est_optimal_ts);
        yhat_suboptimal = predict(m_suboptimal,Proy_est_suboptimal_ts);
        yhat_optimal_t = predict(m_optimal,Proy_est_optimal_t_ts);
        yhat_suboptimal_t = predict(m_suboptimal,Proy_est_suboptimal_t_ts);
 
        yhat_all = predict(mm, vts);
       % yhat_pca = glmval(m_pca, Proy_pca_ts, 'identity');
    elseif strcmpi(model,'logit'),
        
        %clasificando con logistica
        m_optimal = glmfit(Proy_est_optimal_tr,ytr,'binomial');
        m_suboptimal = glmfit(Proy_est_suboptimal_tr,ytr,'binomial');
        m_optimal_t = glmfit(Proy_est_optimal_t_tr,ytr,'binomial');
        m_suboptimal_t = glmfit(Proy_est_suboptimal_t_tr,ytr,'binomial');
        
        %clasificando con naive bayes
        % m_optimal = fitcnb(Proy_est_optimal_tr,ytr,'DistributionNames','kernel');
        % m_suboptimal = fitcnb(Proy_est_suboptimal_tr,ytr,'DistributionNames','kernel') ;
        %m_optimal_t = fitcnb(Proy_est_optimal_t_tr,ytr,'DistributionNames','kernel');
      %  m_suboptimal_t = fitcnb(Proy_est_suboptimal_t_tr,ytr,'DistributionNames','kernel') ;
       
        mm = glmfit(vtr, ytr, 'binomial');
       % m_pca = glmfit(Proy_pca_tr, ytr, 'binomial');
       %prediccion con log?stica
        yhat_optimal(k) = glmval(m_optimal,Proy_est_optimal_ts,'logit')>0.5;
        yhat_suboptimal(k) = glmval(m_suboptimal,Proy_est_suboptimal_ts,'logit')>0.5;
        yhat_optimal_t = glmval(m_optimal,Proy_est_optimal_t_ts,'logit')>0.5;
        yhat_suboptimal_t = glmval(m_suboptimal,Proy_est_suboptimal_t_ts,'logit')>0.5;
        
        %Predicci?n con bayes
        %yhat_optimal = predict(m_optimal,Proy_est_optimal_ts);
        %yhat_suboptimal = predict(m_suboptimal,Proy_est_suboptimal_ts);
        %yhat_optimal_t = predict(m_optimal,Proy_est_optimal_t_ts);
        %yhat_suboptimal_t = predict(m_suboptimal,Proy_est_suboptimal_t_ts);
        
        yhat_all(k) = glmval(mm,vts, 'logit')>0.5;
        %yhat_pca = glmval(m_pca, Proy_pca_ts, 'logit')>0.5;
    end
    
%     err1_optimal(k) = nts-sum(yhat_optimal==yts & yts==1);
%     err2_optimal(k) = nts-sum(yhat_optimal==yts & yts==0);
%     
%     err1_suboptimal(k) = nts-sum(yhat_suboptimal==yts & yts==1);
%     err2_suboptimal(k) = nts-sum(yhat_suboptimal==yts & yts==0);
%     
%     err1_all(k) = nts-sum(yhat_all==yts & yts==1);
%     err2_all(k) = nts-sum(yhat_all==yts & yts==0);

%     err1_nlpca(k) = nts-sum(yhat_nlpca==yts & yts==1);
%     err2_nlpca(k) = nts-sum(yhat_nlpca==yts & yts==0);
%             
%     err_optimal(k) = 1-mean(yhat_optimal==yts);
%     err_suboptimal(k) = 1-mean(yhat_suboptimal==yts);
%     err_optimal_t(k) = 1-mean(yhat_optimal_t==yts);
%     err_suboptimal_t(k) = 1-mean(yhat_suboptimal_t==yts);
%     err_all(k) = 1-mean(yhat_all==yts);
%     err_pca(k) = 1-mean(yhat_pca==yts);


%     
  %  t_optimal_aux(:,:,k) = crosstab(yts, yhat_optimal);
  %  t_suboptimal_aux(:,:,k) = crosstab(yts, yhat_suboptimal);
  %  t_all_aux(:,:,k) = crosstab(yts, yhat_all);
  
%     er_optimal = yts-yhat_optimal;
%     predictions_optimal(k)= er_optimal'*er_optimal/length(yts);
%     pp(k) = sum(er_optimal==0 & yts==1);
%     npnp(k) =sum(er_optimal==0 & yts==0);
%     aciertos_optimal(k) = pp(k) + npnp(k);
%     pnp(k) = sum(er_optimal==1 & yhat_optimal==0);
%     npp(k) =sum(er_optimal ==-1 & yhat_optimal==1);
%     noaciertos_optimal (k) = pnp(k) + npp(k);
%     
%     er_suboptimal = yts-yhat_suboptimal;
%     predictions_suboptimal(k)= er_suboptimal'*er_suboptimal/length(yts);
%     pp(k) = sum(er_suboptimal==0 & yts==1);
%     npnp(k) =sum(er_suboptimal==0 & yts==0);
%     aciertos_suboptimal(k) = pp(k) + npnp(k);
%     pnp(k) = sum(er_suboptimal==1 & yhat_suboptimal==0);
%     npp(k) =sum(er_suboptimal ==-1 & yhat_suboptimal==1);
%     noaciertos_suboptimal (k) = pnp(k) + npp(k);
%     
%     er_all = yts - yhatall;
%     predictions_all(k) = er_all'*er_all/length(yts);
%    pp(k) = sum(er_all==0 & yts==1);
%     npnp(k) =sum(er_all==0 & yts==0);
%     aciertos_all(k) = pp(k) + npnp(k);
%     pnp(k) = sum(er_all==1 & yhatall==0);
%     npp(k) =sum(er_all ==-1 & yhatall==1);
%     noaciertos_all (k) = pnp(k) + npp(k);
%     
%     
%     er_nlpca = yts - yhat_nlpca;
%     predictions_nlpca(k) = er_nlpca'*er_nlpca/length(yts);
%     pp(k) = sum(er_nlpca==0 & yts==1);
%     npnp(k) =sum(er_nlpca==0 & yts==0);
%     aciertos_nlpca(k) = pp(k) + npnp(k);
%     pnp(k) = sum(er_nlpca==1 & yhat_nlpca==0);
%     npp(k) =sum(er_nlpca ==-1 & yhat_nlpca==1);
%     noaciertos_nlpca(k) = pnp(k) + npp(k);
    
end

% misclass1_optimal = mean(err1_optimal);
% misclass2_optimal = mean(err2_optimal);
% 
% misclass1_suboptimal = mean(err1_suboptimal);
% misclass2_suboptimal = mean(err2_suboptimal);
% 
% misclass1_all = mean(err1_all);
% misclass1_all = mean(err2_all);
% 
% misclass1_nlpca = mean(err1_nlpca);
% misclass2_nlpca = mean(err2_nlpca);

misclass_optimal = mean(err_optimal);
misclass_suboptimal = mean(err_suboptimal);

misclass_optimal_t = mean(err_optimal_t);
misclass_suboptimal_t = mean(err_suboptimal_t);
misclass_all = mean(err_all);
%misclass_pca = mean(err_pca);

%   
tabla_optimal = crosstab(Y, yhat_optimal);
tabla_suboptimal = crosstab(Y, yhat_suboptimal);
tabla_all = crosstab(Y, yhat_all);

%  tabla_optimal = nfold*mean(t_optimal_aux, 3);
%  tabla_suboptimal = nfold*mean(t_suboptimal_aux, 3);
%  tabla_all = nfold*mean(t_all_aux, 3);
%     
% mse_optimal = mean(predictions_optimal);
% sd_optimal = std(predictions_optimal);
% aciertosOP =mean(aciertos_optimal);
% NOaciertosOP =mean(noaciertos_optimal);
% 
% %disp('Estimated prediction error from CV procedure is:')
% %disp(mse_optimal)
% 
% 
% mse_suboptimal = mean(predictions_suboptimal);
% sd_suboptimal = std(predictions_suboptimal);
% aciertosSubOP =mean(aciertos_suboptimal);
% NOaciertosSubOP =mean(noaciertos_suboptimal);
% 
% mse_all = mean(predictions_all);
% sd_all = std(predictions_all);
% aciertosALL =mean(aciertos_all);
% NOaciertosALL =mean(noaciertos_all);
% 
% mse_nlpca = mean(predictions_nlpca);
% sd_nlpca = std(predictions_nlpca);
% aciertosNLPCA =mean(aciertos_nlpca);
% NOaciertosNLPCA =mean(noaciertos_nlpca);

%disp('Estimated prediction error from CV procedure is:')
%disp(mse_suboptimal)
