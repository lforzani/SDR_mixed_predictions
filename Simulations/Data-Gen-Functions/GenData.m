function [H, X] = GenData(Gamma,beta,Delta,A,fycent,clase)
% Gamma debe ser un arreglo de qxqxr
% clase = 'mixtos', 'binarios', 'normales'
%GenDataBinariasContinuas_d1(n, q, p, r,Gamma,fycent,Delta,A,beta)

%Generaci?n de datos para el caso de predictores binarios. Modelo para la regreci?n inversa (ising model) con
%H|Y = Gamma1 fy[1]+Gamma2 fy[2]

 %X|H,Y=Afy + beta (Hcent) , varianza Delta.

% muX = zeros(p,1);
% muH = zeros(q,1);
    
 
if (and(nargin == 5, strcmp(clase, 'binarios')))
    n = size(fycent,1);
    p = size(Delta,1);

    H = GenDataBinariasOPv2(n, fycent, Gamma);
    Hcent = H - centering(H);
    for i = 1:n
        X(i,:) = mvnrnd(beta*Hcent(i,:)', Delta);
    end 
    
 elseif (and(nargin == 4, strcmp(clase, 'continuos')))
    n = size(fycent,1);
    p = size(Delta,1);
     
    X = zeros(n,p);
    for i =1:n
      X(i,:) = mvnrnd(A*fycent(i,:)', Delta);
    end
    
elseif (strcmp(clase, 'mixtos'))
    n = size(fycent,1);
    p = size(Delta,1);
     
    H = GenDataBinariasOPv2(n, fycent, Gamma);
    Hcent = centering(H);
    X = zeros(n,p);
    mu1 = zeros(n,p);
    mu2 = zeros(n,p);
    mutotal = zeros(n,p);
    for i = 1:n
      mu1(i,:) = A*fycent(i,:)';
      mu2(i,:) = beta*Hcent(i,:)';
    %   mu2(i,:) = beta*H(i,:)';
      mutotal(i,:) = mu1(i,:)+mu2(i,:);
      X(i,:) = mvnrnd(mutotal(i,:), Delta);
    end
end



 
