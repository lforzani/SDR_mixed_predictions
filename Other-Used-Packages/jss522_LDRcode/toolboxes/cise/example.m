global p;
global d;
global n;
global M;
global N;
nlices = 6;

data = load('D:\boston.txt');

[n1,m] = size(data);

crim = data(1:n1,2);

k = 0;
for i=1:n1
    
    if crim(i) >3.2
        k=k+1;
        ind(k)=i;
    end
    
end

data(ind,:)=[];

[n,m]=size(data);

p = m-1;
d=2;

y = data(:,1);
x = data(:,2:14);
%xr = data(:,2:14);
%xo = xr - ones(n,1)*mean(xr);
%x = xr*diag(1./std(xo));
%[M,N] = SIR(y,x,nlices);
[M,N] = PFC(y,x);

for i=1:101
    
    f = biccis(i/100-0.01);
    stop(i) = f;

end

[f,ind] = min(stop);

%la = fminsearch(@infocis,ind/100-0.01);

lap = ind/100-0.01;


for k = 1:19
    f= biccis(lap-0.01+k/1000);
    dstop(k) = f;
end

[f,dind] = min(dstop);


dlap = lap-0.01+dind/1000;


[fp, betap, stp] = biccis(dlap);


