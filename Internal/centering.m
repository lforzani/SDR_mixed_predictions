function xc = centering(x)
[n,p]=size(x);
xc = x-repmat(mean(x),n,1);