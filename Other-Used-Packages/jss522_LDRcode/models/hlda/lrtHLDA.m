function [Wmin,d,f] = lrtHLDA(Yaux,X,morph,parameters)

%----checking type of response ......................
if strcmpi(morph,'disc'),
    Y = mapdata(Yaux);
    parameters.nslices = max(Y);
end        

% ---- main process............................................................
h = parameters.nslices;
[n,p] = size(X);
data_parameters = setdatapars_v2(Y,X,h);
sigmag = data_parameters.sigmag; 
%--- get handle to objective function, derivative, dof ......................
Fhandle = F(@F4hlda,data_parameters);
dFhandle = dF(@dF4hlda,data_parameters);
dof = @(u) (h*u*(u+1)/2 + (p-u)*(p-u+1)/2 + (h-1)*u + u*(p-u) + p);

f0=0;

ic_choose = icloopHLDA(Fhandle,dFhandle,dof,Y,X,data_parameters,parameters);
d = ic_choose('lrt');
if d>0,
    [Wmin,f] = hlda(Y,X,d,morph,parameters);
else
    f = f0;
    Wmin = zeros(p);
end


