function out = getBin(datos,p)
out = datos;
out.Xtrain = out.Xtrain(:,p+1:end);
out.Xtest = out.Xtest(:,p+1:end);