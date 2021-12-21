%%%This is the script to reproduce figure 2d from the paper by
%%% R. D. Cook and L. Forzani: "Likelihood-based Sufficient Dimension
%%% Reduction". To appear in JASA
%
% BRIEF DESCRIPTION
% The script compares LAD and F2M methods SIR, SAVE and DR by computing the angle 
% between them and the known central subspace. The regression model for the response 
% is Y = X1/4 + a*X1^2/10 + 3*err/5. Figure shows the average angle for different 
% values of parameter 'a' in the regression model. See the paper for details.
% =========================================================================

clear all; 
%setpaths;
nrows = 500;
ncols = 20;
nrep = 20;
amax = 15;

% figure 2d
h = 10;
u = 2;
alp = zeros(ncols,2);
alp(1:3,1) = [1 1 1];
alp(1,2)= 1;
alp(ncols-1,2)=1;
alp(ncols,2)=3;

s0 = zeros(ncols,1);
s0([1 2 3 ncols-1 ncols])=1;
% angulos = zeros(nrep,amax,5);

for a=1:amax
    a
  disp(strcat('a =',int2str(a)));
  for j=1:nrep
    X=normrnd(0,1,nrows,ncols);
 
    yr=.4*a*(X*alp(:,1)).^2 + 3* sin(X*alp(:,2)/4) ;
    y=normrnd(yr,0.2^2);

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
    beta
    angulos(j,a,6)=subspace(beta,alp)*180/pi;
    r1(j,a) = sum(st & s0)/sum(s0);
    r2(j,a) = sum(~st & ~s0)/sum(~s0);
    r3(j,a) = (r1(j,a)*sum(s0) + r2(j,a)*sum(~s0))/length(s0);
  end
end
 
meanang = mean(angulos,1);
save(['pen_jasa2d_' num2str(randint(1,1,[0 5000])) '.mat']);
plot(squeeze(meanang));
%label
title('Y= X_1/4 + aX_2^2/10 + 3\epsilon/5');
xlabel('a');
ylabel('ANGLE');
legend('LAD','SIR','SAVE','DR','AIDA','penAIDA','Location','Best');
