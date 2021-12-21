function [alphan,fn,Wn] = newlad(Yaux,X,dim ,morph,parameters)
% [Wn,fn,fp] = new(Y,X,u,morph
%=====================================================
u = dim(1); d=dim(2); p=cols(X);
%Wlad = load('Wlad.txt');
%----checking type of response and slicing if needed.......................
if strcmpi(morph,'disc'),
    Y = mapdata(Yaux);
    parameters.nslices = max(Y);
else % morph = 'cont'
    if parameters.nslices==0,
        disp('for continuous responses, a number of slices should be given. Five slices will be used');
        parameters.nslices = 5;
    end
    Y = slices(Yaux,parameters.nslices);
end        
%--- get sample statistics ................................................
data_parameters = setdatapars(Y,X,parameters.nslices);
parameters.initvalue = [];

% encuentro estimador inicial para PHI usando HLDA



%--- optimization .........................................................
p = cols(X); aux = eye(u); Wn=aux(1:u,1:d);
% if d == p,
%     disp('WARNING: the subspace you are looking for has the same dimension as the original feature space')
%     fp = Fhandle(Wn);
%     fn = fp;
% else
    count_max = 3;
    stopcount = 0;
    tol = 1e-5; fold = 1e8;
    maxiter = 10; iter = 0;
    [PHIn,ffold] = hlda(Yaux,X,u,'disc',parameters);
    [An,fold] = lad(Yaux,X*PHIn,u,'disc',parameters);
%    An = PHIn'*alphan;
    data_parameters.A = An;

    %--- get handle to objective function and derivative ......................
    Fnewhandle = F(@F4newlad,data_parameters);
    dFnewhandle = dF(@dF4newlad,Fnewhandle);
    dFnewhandle_analitica = dF(@dF4newlad_analitica,data_parameters);
    
    PHIhlda = PHIn; 
    %Wlad = alphan;
    angPHI = []; 
    angLAD = [];
    
%    while (stopcount < count_max)&&(iter<maxiter),
     for iter = 1:5,
%        iter = iter+1;
%        angLAD = [angLAD subspace(Wlad,alphan)*180/pi]
        
        % actualizo PHI
        if isempty(parameters.initvalue)||ischar(parameters.initvalue)
            guess = get_initial_estimate(Y,X,u,data_parameters,parameters);
            PHIo = guess(Fnewhandle); 
            PHIo = orth(PHIo); 
%             PHIo = PHI;
        end
             
        if ~isempty(parameters.sg),
            [fn PHIn] = sg_min(Fnewhandle,dFnewhandle_analitica,PHIo,parameters.sg{:});
        else
            [fn PHIn] = sg_min(Fnewhandle,dFnewhandle_analitica,PHIo,'prcg','euclidean',{1:u},'quiet');
        end
        %angPHI = [angPHI subspace(PHIn,PHIhlda)*180/pi]
        
        % chequea si la F a la salida de sg_min es menor que antes...
        if (fn > Fnewhandle(PHIo)),
            error('SG_MIN NO ESTA MINIMIZANDO!!!');
        end
        
        PHIX = X*PHIn;
        
        % actualizo A;
        [An,ffn,ddd] = lad(Y,PHIX,d,'disc',parameters);
        
 %       if abs(ffn-fold)<tol,
  %         stopcount = stopcount+1;
   %     end
        fold = ffn;
        alphan = PHIn*An;
        data_parameters.A = An;
        [Fnewhandle, dFnewhandle] = update_handles(data_parameters);

    end
    alphan = PHIn*An;
    Wn = PHIn;
% end

