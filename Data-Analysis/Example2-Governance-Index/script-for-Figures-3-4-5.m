%This script is to run the Figures in the example of Governance index.
%

clear all;

%set seed


%Open-Read data set
datos = readtable('~/Data-Analysis/Example2-Governance-Index/GovernanceDataExtended2.csv', 'Delimiter', ';', 'HeaderLines',1);


% Specify the number of continuous (p) , binary predictors (q) y r 
p = 6;
q=11;
r=1;
X = table2array(datos(:,2:7));
H = table2array(datos(:,10:(10+(q-1))));
Y = table2array(datos(:,1));

%% For Figures 4, 5 and 6
%First we obtain CG index for overall sample using PCA, PFC and SDRoptimal
%PCA

[coef, scores] = pca(X);
CGindex_PCA = scores(:,1); 

%PFC
Fy = get_fy(Y,r);
[Proy_est_PFC,redu_est_PFC] = ldr(Y,X,'pfc','cont',1,'Fy',Fy);
%CGindex_PFC=Proy_est_PFC;
CGindex_PFC=X*(-1)*redu_est_PFC;

%Optimal and SubOptimal

[red_est_optimal] = rr4mix_optimal(X,H,Y,1,r,'auto');

 [red_est_suboptimal1,red_est_suboptimal2] = rr4mix_suboptimal(X,H,Y,1,1,r,'auto');
       

k=size(H,2)*(size(H,2)+1)/2;
n=size(H,1);
    Hraro = zeros(n,k);
    for j = 1:n
       aux = H(j,:)'*H(j,:);
      Hraro(j,:) = [diag(aux);aux(find(tril(ones(size(H,2),size(H,2)),-1)))];
    end
    T_optimal = [X,Hraro];

    CGindex_SDROptimal  = T_optimal*red_est_optimal;

CGIndexes=[CGindex_PCA CGindex_PFC CGindex_SDROptimal];

%cd '~\Data_Analysis\Example2-Governance-Index';

%csvwrite('CGIndexes.csv',CGIndexes);

% Read PCAmix from R
PCAmix =readtable('~\Data-Analysis\Example2-Governance-Index\PCAmixIndex.csv', 'Delimiter', ',');
CGindex_PCAmix= table2array(PCAmix(:,2));

%Plot fit CGindex by PCA and PFC
rPCA= ksr(CGindex_PCA,Y);
rPFC= ksr(CGindex_PFC,Y);

lfitPCA = fitlm(CGindex_PCA,Y);
yhatpca= predict(lfitPCA, CGindex_PCA);

lfitPFC = fitlm(CGindex_PFC,Y);
yhatpfc= predict(lfitPFC, CGindex_PFC);

graf1=gscatter(CGindex_PCA,Y)
hold on;
graf2=plot(rPCA.x,rPCA.f,'b',CGindex_PCA, yhatpca,'k');
 ylabel('ln GDP');
    xlabel('CGindex-PCA');
legend(graf2,{'Kernel', 'Linear'},'NumColumns', 2);

graf3=gscatter(CGindex_PFC,Y)
hold on;
graf4=plot(rPFC.x,rPFC.f,'b',CGindex_PFC, yhatpfc,'k');
 ylabel('ln GDP');
    xlabel('CGindex-PFC');
legend(graf4,{'Kernel', 'Linear'},'NumColumns', 2);

%Countries as factors
factor=table2array(datos(:,8));

%PCAmix
lfitPCAmix = fitlm(    CGindex_PCAmix,Y);
yhatpcamix= predict(lfitPCAmix,     CGindex_PCAmix);
rPCAmix= ksr(CGindex_PCAmix,Y);
graf5=gscatter(CGindex_PCAmix,Y,factor)
hold on;
graf6=plot(rPCAmix.x,rPCAmix.f,'b',CGindex_PCAmix, yhatpcamix,'k');
 ylabel('ln GDP');
    xlabel('CG-PCAmix');
legend(graf5, graf6,{'Kernel', 'Linear'});


%SDR Optimal
lfitPFCmix = fitlm(CGindex_SDROptimal,Y);
yhatpfcmix= predict(lfitPFCmix,     CGindex_SDROptimal);
rPFCmix= ksr(CGindex_SDROptimal,Y);

graf6=gscatter(CGindex_SDROptimal,Y,factor)
hold on;
graf7=plot(rPFCmix.x,rPFCmix.f,'b',CGindex_SDROptimal, yhatpfcmix,'k');
 ylabel('ln GDP');
    xlabel('CG-SDROptimal');
legend(graf6, graf7,{'Kernel', 'Linear'});

