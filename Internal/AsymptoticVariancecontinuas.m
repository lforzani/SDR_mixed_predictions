function  Vrcl = AsymptoticVariancecontinuas(Deltahat,Ahat,fycent,muX);
    % RODRIGO: ver los argumentos que tenemos que pedir
    % AsymptoticVariancecontinuas(Deltahat,Ahat,fycent,muX)
    
    % Matriz de varianza asintotica para el test de dimensi?n

    n = size(fycent,1);
    p = size(Ahat,1);
    r = size(fycent,2);

    if (nargin == 3)
        muX = zeros(p,1);
    end

    Ip = eye(p);
    Ir = eye(r);
    Ipr = eye(p*r);
    Ipp = eye(p*(p+1)/2);
    
    
    invDelta = inv(Deltahat);
    
    B = invDelta*[muX, Ahat]; 
    
    K = inv(vecperm(p,r+1));
    uno = ones(n,1);
    Fy = [uno, fycent];
    Sigma = Fy'*Fy/n;
    invSigma = inv(Sigma);
    V = kron(invSigma, invDelta)+ kron(B'*Deltahat*B,invDelta) + kron(B',B)*K;


    
   
trC = [zeros(r,1) , Ir];
Vrcl = kron(trC,Ip)*V*kron(trC',Ip);
end



