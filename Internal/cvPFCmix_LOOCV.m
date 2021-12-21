function [MSE_optimal, sd_optimal, MSE_suboptimal, sd_suboptimal, MSE_all, sd_all, MSE_PFC, sd_PFC, ProyPFC_tr,ProyPFC_ts, ProyOptimal_tr, ProyOptimal_ts,ProySuboptimal_tr1, ProySuboptimal_ts1,ProySuboptimal_tr2, ProySuboptimal_ts2] = cvPFCmix_LOOCV(Y,X,H,type,model)

rng(1000)

disp('== starting CV procedure for PFCMix (this can take several minutes)===')

% this iis an SES index so we dont need to choose the dimension
dim = 1;
dim1 = 1;
dim2 = 1;    

% thesholds for cutoff
trh = 0.1:0.05:0.9;

n = size(Y,1);
p = size(X,2);
q = size(H,2);
kk = q*(q+1)/2;
r = type; 

error_optimal = zeros(1,n);
%lambdas_optimal = zeros(1,n);
error_suboptimal = zeros(1,n); 
error_all = zeros(1,n);
error_PFC = zeros(1,n);
ProyPFC_tr = zeros(n-1,n);
ProyPFC_ts = zeros(1,n);
ProyOptimal_tr = zeros(n-1,n);
ProyOptimal_ts = zeros(1,n);
ProySuboptimal_tr1 = zeros(n-1,n);
ProySuboptimal_ts1 = zeros(1,n);
ProySuboptimal_tr2 = zeros(n-1,n);
ProySuboptimal_ts2 = zeros(1,n);

ProyOptimal0_tr = zeros(n-1,n);
ProyOptimal0_ts = zeros(1,n);
ProySuboptimal0_tr1 = zeros(n-1,n);
ProySuboptimal0_ts1 = zeros(1,n);
ProySuboptimal0_tr2 = zeros(n-1,n);
ProySuboptimal0_ts2 = zeros(1,n);


eroptimal =zeros(1,n);
eroptimal0 =zeros(1,n);
ersuboptimal=zeros(1,n);
ersuboptimal0=zeros(1,n);


for k = 1:n
    rng(1000)
    k
    xts = X(k,:); yts = Y(k,:); hts = H(k,:);
    xtr = X; xtr(k,:) = []; ytr = Y; ytr(k)=[]; htr = H; htr(k,:)=[];
    nts = size(xts,1);
    ntr = size(xtr,1);
      
    %estimate the reduction with optimal and subopitmal PFCmix
    [red_est_optimal] = rr4mix_optimal(xtr,htr,ytr,dim,type,'auto');
    [red_est_suboptimal1,red_est_suboptimal2] = rr4mix_suboptimal(xtr,htr,ytr,dim1,dim2,type,'auto');
       
    %without regularization
    red_est_optimal0 = rr4mix_optimal(xtr,htr,ytr,dim,type,0);
    [red_est_suboptimal01,red_est_suboptimal02] = rr4mix_suboptimal(xtr,htr,ytr,dim1,dim2,type,0);
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
    %Proy_est_PFC_tr = Proy_est_PFC_tr;
    Proy_est_optimal_tr = T_optimal_tr*red_est_optimal;
    Proy_est_optimal_ts = T_optimal_ts*red_est_optimal;
    Proy_est_suboptimal_tr1 = T_suboptimal_tr1*red_est_suboptimal1; 
    Proy_est_suboptimal_ts1 = T_suboptimal_ts1*red_est_suboptimal1;
          
    
    Proy_est_optimal0_tr = T_optimal_tr*red_est_optimal0;
    Proy_est_optimal0_ts = T_optimal_ts*red_est_optimal0;
    Proy_est_suboptimal0_tr1 = T_suboptimal_tr1*red_est_suboptimal01; 
    Proy_est_suboptimal0_ts1 = T_suboptimal_ts1*red_est_suboptimal01;
    
    if(size(red_est_suboptimal2,2) == 0)
        Proy_est_suboptimal_tr = Proy_est_suboptimal_tr1;
        Proy_est_suboptimal_ts = Proy_est_suboptimal_ts1;
    else     
       Proy_est_suboptimal_tr2 = T_suboptimal_tr2*red_est_suboptimal2;
        Proy_est_suboptimal_ts2 = T_suboptimal_ts2*red_est_suboptimal2; 

        Proy_est_suboptimal_tr = [Proy_est_suboptimal_tr1,Proy_est_suboptimal_tr2];
        Proy_est_suboptimal_ts = [Proy_est_suboptimal_ts1,Proy_est_suboptimal_ts2];
    end
    
        
    if(size(red_est_suboptimal02,2) == 0)
        Proy_est_suboptimal0_tr = Proy_est_suboptimal0_tr1;
        Proy_est_suboptimal0_ts = Proy_est_suboptimal0_ts1;
    else     
       Proy_est_suboptimal0_tr2 = T_suboptimal_tr2*red_est_suboptimal02;
        Proy_est_suboptimal0_ts2 = T_suboptimal_ts2*red_est_suboptimal02; 

        Proy_est_suboptimal0_tr = [Proy_est_suboptimal0_tr1,Proy_est_suboptimal0_tr2];
        Proy_est_suboptimal0_ts = [Proy_est_suboptimal0_ts1,Proy_est_suboptimal0_ts2];
    end
    
    
    %Save the proyections to use in R
    %ProyOptimal_tr(:,k) = Proy_est_optimal_tr;
    %ProyOptimal_ts(:,k) = Proy_est_optimal_ts;
    % ProySuboptimal_tr1(:,k) = Proy_est_suboptimal_tr1;
    % ProySuboptimal_ts1(:,k) = Proy_est_suboptimal_ts1;
    % ProySuboptimal_tr2(:,k) = Proy_est_suboptimal_tr2;
    % ProySuboptimal_ts2(:,k) = Proy_est_suboptimal_ts2;

    
   if strcmpi(model,'linear')
       rng(1000)
        %Here we include PFC fit for only continuous variables
        %fit the model with training
        % estimate the reduction with PFC
        %if ~isstr(r),
         %   Fy = get_fy(Y,r);
        %else
         %   Fy = get_fyZ(Y);
        %end
        Fy = get_fy(ytr,r);
        [Proy_est_PFC_tr,redu_est_PFC] = ldr(ytr,xtr,'pfc','cont',r,'Fy',Fy);
        Proy_est_PFC_ts = xts*redu_est_PFC;
        ProyPFC_tr(:,k) = Proy_est_PFC_tr;
        ProyPFC_ts(:,k) = Proy_est_PFC_ts;
       
        m_optimal =  fitlm(Proy_est_optimal_tr,ytr);
        m_suboptimal =  fitlm(Proy_est_suboptimal_tr,ytr);
        m_optimal0 =  fitlm(Proy_est_optimal0_tr,ytr);
        m_suboptimal0 =  fitlm(Proy_est_suboptimal0_tr,ytr);
        m_PFC = fitlm(Proy_est_PFC_tr,ytr);
        mm = fitlm(vtr, ytr);
        
        % predict testing   
        yhat_optimal = predict(m_optimal,Proy_est_optimal_ts); 
        yhat_suboptimal = predict(m_suboptimal,Proy_est_suboptimal_ts);
        yhat_optimal0 = predict(m_optimal0,Proy_est_optimal0_ts); 
        yhat_suboptimal0 = predict(m_suboptimal0,Proy_est_suboptimal0_ts);
        
        yhat_PFC = predict(m_PFC,Proy_est_PFC_ts);
        yhat_all = predict(mm, vts);
        
        eroptimal(k) = mean((yhat_optimal-yts).^2);
        eroptimal0(k) = mean((yhat_optimal0-yts).^2);

        ersuboptimal(k) = mean((yhat_suboptimal-yts).^2);
        ersuboptimal0(k) = mean((yhat_suboptimal0-yts).^2);
        
        error_optimal(k) =min(eroptimal(k),eroptimal0(k));
        error_suboptimal(k)=min(ersuboptimal(k), ersuboptimal0(k));
        
        %if error_optimal(k) == eroptimal(k),
          %  lambdas_optimal(k) =lambdaopt;
        %else
          %  lambdas_optimal(k) = 0;
        %end
        
        error_PFC(k) = mean((yhat_PFC-yts).^2);
        error_all(k) = mean((yhat_all-yts).^2);
        
        if error_optimal(k) == eroptimal(k),
             ProyOptimal_tr(:,k) = Proy_est_optimal_tr;
             ProyOptimal_ts(:,k) = Proy_est_optimal_ts;
        else
             ProyOptimal_tr(:,k) = Proy_est_optimal0_tr;
             ProyOptimal_ts(:,k) = Proy_est_optimal0_ts;
        end
        
        if  error_optimal(k) == eroptimal(k),
               ProySuboptimal_tr1(:,k) = Proy_est_suboptimal_tr1;
               ProySuboptimal_ts1(:,k) = Proy_est_suboptimal_ts1;
               ProySuboptimal_tr2(:,k) = Proy_est_suboptimal_tr2;
               ProySuboptimal_ts2(:,k) = Proy_est_suboptimal_ts2;
        else
             ProySuboptimal_tr1(:,k) = Proy_est_suboptimal0_tr1;
               ProySuboptimal_ts1(:,k) = Proy_est_suboptimal0_ts1;
               ProySuboptimal_tr2(:,k) = Proy_est_suboptimal0_tr2;
               ProySuboptimal_ts2(:,k) = Proy_est_suboptimal0_ts2;
        end
        
   elseif strcmpi(model,'logit')
       rng(1000)
        % fit the model with training
        m_optimal = glmfit(Proy_est_optimal_tr,ytr,'binomial');
        m_suboptimal = glmfit(Proy_est_suboptimal_tr,ytr,'binomial');
        m_optimal0 = glmfit(Proy_est_optimal0_tr,ytr,'binomial');
        m_suboptimal0 = glmfit(Proy_est_suboptimal0_tr,ytr,'binomial');

        mm = glmfit(vtr, ytr,'binomial');

        % choose the optimal cut off 
        c_optimal = cutoff_optimo_minMSE(m_optimal,ytr,Proy_est_optimal_tr,trh);
        c_suboptimal = cutoff_optimo_minMSE(m_suboptimal,ytr,Proy_est_suboptimal_tr,trh);
        
        c_optimal0 = cutoff_optimo_minMSE(m_optimal0,ytr,Proy_est_optimal0_tr,trh);
        c_suboptimal0 = cutoff_optimo_minMSE(m_suboptimal0,ytr,Proy_est_suboptimal0_tr,trh);
        
        
        c_all = cutoff_optimo_minMSE(mm,ytr,vtr,trh);

        % predict testing with the optimal cut off    
        yhat_optimal = glmval(m_optimal,Proy_est_optimal_ts,'logit')>c_optimal;
        yhat_suboptimal = glmval(m_suboptimal,Proy_est_suboptimal_ts,'logit')>c_suboptimal;
        
        yhat_optimal0 = glmval(m_optimal0,Proy_est_optimal0_ts,'logit')>c_optimal0;
        yhat_suboptimal0 = glmval(m_suboptimal0,Proy_est_suboptimal0_ts,'logit')>c_suboptimal0;
        
        yhat_all = glmval(mm,vts,'logit')>c_all;
        
        
        % compute errors   
        
        eroptimal(k) = mean((yhat_optimal-yts).^2);
        eroptimal0(k) = mean((yhat_optimal0-yts).^2);

        ersuboptimal(k) = mean((yhat_suboptimal-yts).^2);
        ersuboptimal0(k) = mean((yhat_suboptimal0-yts).^2);
        
        error_optimal(k) =min(eroptimal(k),eroptimal0(k));
        error_suboptimal(k)=min(ersuboptimal(k), ersuboptimal0(k));
        
        error_all(k) = mean((yhat_all-yts)^2);
   %      if error_optimal(k) == eroptimal(k),
   %         lambdas_optimal(k) =lambdaopt;
  %%      else
          %  lambdas_optimal(k) = 0;
  %      end
        
      %  error_PFC(k) = mean((yhat_PFC-yts).^2);
        error_all(k) = mean((yhat_all-yts).^2);
       
        
        
        if error_optimal(k) == eroptimal(k),
             ProyOptimal_tr(:,k) = Proy_est_optimal_tr;
             ProyOptimal_ts(:,k) = Proy_est_optimal_ts;
        else
             ProyOptimal_tr(:,k) = Proy_est_optimal0_tr;
             ProyOptimal_ts(:,k) = Proy_est_optimal0_ts;
        end
        
        if  error_optimal(k) == eroptimal(k),
               ProySuboptimal_tr1(:,k) = Proy_est_suboptimal_tr1;
               ProySuboptimal_ts1(:,k) = Proy_est_suboptimal_ts1;
               ProySuboptimal_tr2(:,k) = Proy_est_suboptimal_tr2;
               ProySuboptimal_ts2(:,k) = Proy_est_suboptimal_ts2;
        else
             ProySuboptimal_tr1(:,k) = Proy_est_suboptimal0_tr1;
               ProySuboptimal_ts1(:,k) = Proy_est_suboptimal0_ts1;
               ProySuboptimal_tr2(:,k) = Proy_est_suboptimal0_tr2;
               ProySuboptimal_ts2(:,k) = Proy_est_suboptimal0_ts2;
        end
    end        
end

MSE_optimal = mean(error_optimal);
MSE_suboptimal = mean(error_suboptimal);
MSE_PFC = mean(error_PFC);
MSE_all = mean(error_all);


sd_optimal = std(error_optimal);
sd_suboptimal = std(error_suboptimal);
sd_PFC = std(error_PFC);
sd_all = std(error_all);
