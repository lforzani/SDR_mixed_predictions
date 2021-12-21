function out = getCont(datos,p)
out = datos;
out.Xtrain = out.Xtrain(:,1:p);
out.Xtest = out.Xtest(:,1:p);