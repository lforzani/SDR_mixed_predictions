%%%This is the script to reproduce figure 2b from the paper by
%%% R. D. Cook and L. Forzani: "Likelihood-based Sufficient Dimension
%%% Reduction". To appear in JASA
%
% BRIEF DESCRIPTION
% The script compares LAD and F2M methods SIR, SAVE and DR by computing the angle 
% between them and the known central subspace. The regression model for the response 
% is Y = X1^2/(20*a) + err/10. Figure shows the average angle for different values of 
% parameter 'a' in the regression model. See the paper for details.
% =========================================================================


clear all; 
%setpaths;
nrows = 500;
ncols = 8;
nrep = 50;
amax = 20;

% figure 2b
h = 5;
u = 1;
alp = zeros(ncols,1);
alp(1,1) = 1;


angulos = zeros(nrep,amax,5);

for a=1:amax
  disp(strcat('a =',int2str(a)));
  for j=1:nrep
    X=normrnd(0,1,nrows,ncols);
 
%   figure 2b  
   yr=(X*alp(:,1)).^2/(20*a);
    y=normrnd(yr,0.1^2);
 
    [WX, W]=ldr(y,X,'lad','cont',u,'nslices',h);
    angulos(j,a,1)=subspace(W,alp)*180/pi;

    [WX, W]=SIR(y,X,'cont',u,'nslices',h);
    angulos(j,a,2)=subspace(W,alp)*180/pi;
    
    [WX, W]=SAVE(y,X,'cont',u,'nslices',h);
    angulos(j,a,3)=subspace(W,alp)*180/pi;

    [WX, W]=DR(y,X,'cont',u,'nslices',h);    
    angulos(j,a,4)=subspace(W,alp)*180/pi;
    
     [Wx,W]=aida(y,X,u,'cont');
    angulos(j,a,5)=subspace(W,alp)*180/pi;
    
%     nfolds=10;
%     lambdas = .001:.005:1;
    [f,beta,st]=cise(y,X,u,1,'AIDA');
    angulos(j,a,6)=subspace(beta,alp)*180/pi;
    r1(j,a) = sum(beta & alp)/sum(alp);
    r2(j,a) = sum(~beta & ~alp)/sum(~alp);
    r3(j,a) = (r1(j,a)*sum(alp) + r2(j,a)*sum(~alp))/length(alp);
  end
end

meanang = mean(angulos,1);
save(['pen_jasa2b_' num2str(randint(1,1,[0 5000])) '.mat']);

plot(squeeze(meanang));
%label
title('Y=X_1^2/(20a) + \epsilon/10');
xlabel('a');
ylabel('ANGLE');
legend('LAD','SIR','SAVE','DR','AIDA','penAIDA','Location','Best');
