function [H] = GenDataBinariasOPv2(n, fycent, Gamma)
%Generación de datos para el caso de predictores binarios. Modelo para la regreci?n inversa (ising model) con
	%H|Y = Gamma1 fy[1]
	
% Constantes ---------------------------------------------------------------------------------------
		[q, ~, r] = size(Gamma);

% Generacion de discretas --------------------------------------------------------------------------
		H = zeros(n,q);
	
		for i = 1:n
		
				aux = zeros(q,q);
			
				for t=1:r
                    auxxx = reshape(Gamma(:,:,t),[q,q]);
					aux = aux+auxxx*fycent(i,t); %Equivale a hacer = Gamma1*auu1(s,1)+Gamma2*auu1(s,2)+ ...
				end
			
				aux1 = diag(aux);	
				H(i,:) = binornd( ones(size(aux1)), 1./(1+exp(-aux1)), size(aux1));
			
		end
	
% Eliminar filas repatidas --------------------------------------------------------------------------
		[auu1, ia, ic]  = unique(fycent,'rows'); 
	    r_auu1=size(auu1,1);
	
		for s=1:r_auu1
				nn(s) = sum(ic == ic(ia(s)) ); %contar cantidad de filas repetidas en fycent
		end
	
	  nreps = 1000;  
	
  	for j=1:nreps;
		for s=1:r_auu1 %(r+1);
			
			aux = zeros(q,q);
			for t=1:r
					aux=aux+ reshape(Gamma(:,:,t),[q,q])*auu1(s,t); %Equivale a hacer = Gamma1*auu1(s,1)+Gamma2*auu1(s,2)+ ...
			end	
			
			etass=reshape(aux,[q, q]);
			
			for hh=1:q; 
					Yraro = H(ic == ic(ia(s)), :);
					Yraro(:, hh) = 1;
					etas = Yraro*etass;
				
					ss=etas(:,hh);
					H(ic == ic(ia(s)), hh) = binornd(ones(size(ss)), 1./(1+exp(-ss)),nn(s),1);
			end
	
		end
  	end
	
%	tauhat = levina4ising(H,fycent,0.001);
% Comparar con Gamma1
%	reshape(tauhat(2,:,:),[q q]);
%	auxx2 = ones(q,q);
%	TrueRedBin = [diag(Gamma1); Gamma1(find(tril(auxx2,-1)))];
	
end
