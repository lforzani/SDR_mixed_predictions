%%% This is Figure 2a from R. D. Cook paper: Fisher Lecture: Dimension
%%% Reduction in Regression. Statist. Sci. Volume 22, Number 1 (2007), 1-26
%
% It aims to compare likelihood-based estimate to other methods and to assess the 
% effect of initial candidates on the obtained estimate. In this subfigure,
% the PFC estimate is obtained from the PC candidate set. The average angle
% between the estimate and the true central subspace is taken as a measure
% of performance.
% We reproduce here OLS, SIR and the actual MLE in Grassmann manifold. 
% Therefore this picture it is a little different than the original one.

clear all; close all;
%setpaths; % adds all directories in the LBDR package to MATLAB's path. 
randn('state',10);
rand('state',10);
nrows=500;
sigg=1:.5:4.1;
sigy=1;
sig0=1;
ncols=10; nreps=25; u=1;
ang=zeros(5,length(sigg),nreps);

lambdas = .001:.005:1;
nfolds = 5;

for k=1:length(sigg)
    disp(k);
    sig=sigg(k);
    alp=zeros(ncols,1);
    alp(1)=1;
    alp0=[zeros(ncols-1,1)'
    eye(ncols-1)];
    for j=1:nreps
       Y=zeros(nrows,1);
       X=zeros(nrows,ncols);
       for hh=1:nrows
          Y(hh)=normrnd(0,sigy);
          t1 = normrnd(0,1);
          t0 = mvnrnd(zeros(ncols-1,1),eye(ncols-1));
          X(hh,:)=  (sig*alp*t1)'+(sig0*alp0*t0')'   + Y(hh)*alp';
       end
      [WX,W]=SIR(Y,X,'cont',1,'nslices',8);
      ang(1,k,j)=subspace(W,alp)*180/pi;
      [W]=  regress(Y,X);
      ang(2,k,j)=subspace(W,alp)*180/pi;
      FF = get_fy(Y,1);
      [WX,W] = ldr(Y,X,'EPFC','cont',1,'fy',FF);
      ang(3,k,j)=subspace(W,alp)*180/pi;
      [WX,W] = eigEPFC(Y,X,1,'cont',FF);
      ang(4,k,j)=subspace(W,alp)*180/pi;
      [WX,W] = ldr(Y,X,'EPFC','cont',1,'fy',FF,'initval',W);
      ang(5,k,j)=subspace(W,alp)*180/pi;
      
      [f,W,st] = cise(Y,X,1,1,'EPFC');
      ang(6,k,j)=subspace(W,alp)*180/pi;
      beta=W;
    r1(k,j) = sum(beta & alp)/sum(alp);
    r2(k,j) = sum(~beta & ~alp)/sum(~alp);
    r3(k,j) = (r1(k,j)*sum(alp) + r2(k,j)*sum(~alp))/length(alp);
      
      [WX,W] = ldr(Y,X,'EPFC','cont',1,'fy',FF,'initval',W);
      ang(7,k,j)=subspace(W,alp)*180/pi;

    end
end
  
angmean=zeros(3,length(sigg));

for k=1:length(sigg)
   for p=1:5
      angmean(p,k) = mean(ang(p,k,:) );
   end
end

plot(sigg,squeeze(angmean(1,:)),'-',sigg,squeeze(angmean(2,:)),'-',sigg,squeeze(angmean(3,:)),'-',sigg,squeeze(angmean(4,:)),'-',sigg,squeeze(angmean(5,:)),'-',sigg,squeeze(angmean(6,:)),'-',sigg,squeeze(angmean(7,:)),'-')
% labels
xlabel('\sigma'); ylabel('Angle');
ylim([0 80]);
legend('SIR','OLS','EPFC','eigEPFC','eigINIT','CISEepfc','peninitEPFC');
title('Angle vs. \sigma');
