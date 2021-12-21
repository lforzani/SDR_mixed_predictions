function [Wn,fn,fp] = hlda(Yaux,X,u,morph,parameters)
% [Wn,fn,fp] = lad(Y,X,u,morph,parameters)

%----checking type of response ......................
if strcmpi(morph,'disc'),
    Y = mapdata(Yaux);
    parameters.nslices = max(Y);
else % morph = 'cont'
    Y = Yaux;
    parameters.nslices = length(Y);
end        

%--- get sample statistics ................................................
data_parameters = setdatapars(Y,X,parameters.nslices);

if strcmpi(morph,'cont')
    SIGMAfit = get_fitted_cov(Y,X,parameters.fy);
else
    SIGMAfit = get_average_cov(X,data_parameters);
end
SIGMA = data_parameters.sigmag;
SIGMAres = SIGMA - SIGMAfit;
data_parameters.Afit = SIGMAfit;
data_parameters.B = SIGMAres;


%--- get handle to objective function and derivative ......................
Fhandle = F(@F4hlda,data_parameters);
dFhandle = dF(@dF4hlda,data_parameters);

%--- get initial estimate .................................................
if strcmpi(morph,'cont')
    haux = 5;
    Ysliced = slices(Y,haux);
    aux_datapars = setdatapars(Ysliced,X,haux);
    auxpars = parameters; auxpars.nslices=haux;
else
    Ysliced = Y;
    aux_datapars = data_parameters;
    auxpars = parameters;
end

if isempty(parameters.initvalue)||ischar(parameters.initvalue)
    guess = get_initial_estimate(Ysliced,X,u,aux_datapars,auxpars);
    Wo = guess(Fhandle);
else
    Wo = parameters.initvalue;
end

%--- optimization .........................................................
p = cols(X); Wn = eye(p);
fp = Fhandle(Wn);
if u == p,
    disp('WARNING: the subspace you are looking for has the same dimension as the original feature space')
    fn = fp;
else
    if ~isempty(parameters.sg),
        [fn Wn] = sg_min(Fhandle,dFhandle,Wo,parameters.sg{:});
    else
        [fn Wn] = sg_min(Fhandle,dFhandle,Wo,'prcg','euclidean',{1:u},'quiet');
    end
end

