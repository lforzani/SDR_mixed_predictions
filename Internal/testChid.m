function [d1,d2] = testChid(n,Vrcl,bhat,alpha)

    m = size(bhat,1);
    r = size(bhat,2);
    k = min(m,r);

    d = rank(bhat);

    % Test for d (Weighted chi-square test)

    [U,D,V] = svd(bhat);
    
    for j = 0:(k-1) 

        s = min(rank(Vrcl), (r-j)*(m-j));

        U0 = U(:,(j+1):m);
        R0 = V(:,(j+1):r); 
        K0 = D((j+1):m,(j+1):r);

        Qhat = kron(R0',U0')*Vrcl*kron(R0,U0);

        Lambda1 = n*(vec(K0))'*vec(K0);

        %-------------------------------
        % permutate
        N = 50000;

        w = svd(Qhat);

        cont = 0;
        for i = 1:N 
              aux = 0;

              for dd=1:rank(Qhat)
                  aux = aux + w(dd)*chi2rnd(1);
              end

              if aux > Lambda1
                  cont = cont + 1;
              end
        end

        pvalor1 = cont/N;

        %-------------------------------

        if (pvalor1 > alpha)
            % disp('Do not reject H0')
            d1 = j;
            break
        else
            d1 = d;
        end 
    end
    
    [U,D,V] = svd(bhat);

    % Chi-square test version 1
    d2=d;
    
    for j = 0:(k-1) %ver esto
        U0 = U(:,(j+1):m);
        R0 = V(:,(j+1):r); 
        K0 = D((j+1):m,(j+1):r);

        s = min(rank(Vrcl), (r-j)*(m-j));

        Qhat = kron(R0',U0')*Vrcl*kron(R0,U0);

        Lambda2 = n*(vec(K0))'*pinv(Qhat)*vec(K0);

        [QQ MM VV]=svd(Qhat);
        S=svd(Qhat);
        S(1:rank(Qhat))=1./S(1:rank(Qhat));
        S((rank(Qhat)+1):size(Qhat,2))=0;
        cc=QQ*diag(S)*VV';
        Lambda2 = n*(vec(K0))'*cc*vec(K0);

        pvalor2 = 1 - chi2cdf(Lambda2,s);

        if (pvalor2 > alpha)
            % disp('Do not reject H0')
            d2 = j;
        break
    end
end



