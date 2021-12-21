function [Wmin,d,fmin] = permLAD(Yaux,X,morph,parameters)
%----checking type of response and slicing if needed.......................
if strcmpi(morph,'disc'),
    Y = mapdata(Yaux);
    parameters.nslices = max(Y);
end
h = parameters.nslices;
alpha = parameters.alpha;
nsample = parameters.nsamples;

[n,p] = size(X);
data_parameters = setdatapars_v2(Y,X,h);

F0handle = F(@F4hlda,data_parameters);
dF0handle = dF(@dF4hlda,data_parameters);
% f0 = @(sigmag) (n*p*(1+log(2*pi))/2 + n*logdet(sigmag)/2);

for u=1:p-1.
    u
    fpo = F0handle(eye(p));
    guess = get_initial_estimate(Y,X,u,data_parameters,parameters);
    Wo = guess(F0handle);
    [fno,Wu] = sg_min(F0handle,dF0handle,Wo,'prcg','euclidean',{1:u},'quiet');
    T = 2*(fno-fpo);

    [Q R] = qr(Wu);
    Wu0 = Q(:,(u+1):p);
    Xnew = X*[Wu Wu0];
    Xrsp = zeros(n,p);
    fns = 1:nsample;
    fps = 1:nsample;
    for i=1:nsample,
        YY = randperm(n); % ver implementación alternativa en denboot.m
        Xrsp(:,1:u) = Xnew(:,1:u);
        Xrsp(:,(u+1):p) = Xnew(YY(:),(u+1):p);
        newpars = setdatapars_v2(Y,Xrsp(:,:),h);
        Fnew = F(@F4hlda,newpars);
        dFnew = dF(@dF4hlda,newpars);
        fps(i) = Fnew(eye(p));
        guessnew = get_initial_estimate(Y,X,u,newpars,parameters);
        Wo = guessnew(Fnew);
        fns(i) = sg_min(Fnew,dFnew,Wo,'prcg','euclidean',{1:u},'quiet');
    end
    fnmfp = 2*(fns - fps);
    prop = sum(fnmfp > T);
    aux = prop/nsample;
    if (aux > (alpha))||(u==p-1),
        d = u;
        [Wmin,fmin] = hlda(Y,X,d,morph,parameters);
        break;
    end
end

 
