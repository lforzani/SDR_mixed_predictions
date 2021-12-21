function cutoff_optimal = cutoff_optimo_minMSE(modelo_ajustado,Y,X,trh)
 
ntrh = size(trh, 2);

error = zeros(ntrh,1);

for i = 1:ntrh
    yhat = glmval(modelo_ajustado,X,'logit')>trh(i);
   
    error(i) = 1-mean(yhat==Y);
end

[~,idx] = min(error);
cutoff_optimal = trh(idx);
