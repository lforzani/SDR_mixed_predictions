function [b,fval] = fzero_vDT(fun,x,options,Y,dr)
% Descripcion minima:
% Inputs:
% - fun: handle a la función a la que se quiere encontrar una raiz. 
% Dados Y y dr, fun es una función solamente de x y buscamos la raiz 
% en el intervalo [x(1); x(2)]
% - options no lo usamos. Lo puse solamente por compatibilidad
% Outputs
% - b: raíz de fun en el intervalo [x(1),x(2)]
% - fval = fun(b)


%%
% Initialization
fcount = 0;
iter = 0;
intervaliter = 0;
exitflag = 1;
procedure = ' ';

tol = eps;


% Interval input
if (length(x) == 2) 
    a = x(1); savea=a;
    b = x(2); saveb=b;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    fa = fun(a,Y,dr);
    fb = fun(b,Y,dr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    fcount = fcount + 2;
    savefa = fa; savefb = fb;
    
    if ( fa == 0 )
        b = a;
        fval = fa;
        return
    elseif ( fb == 0)
        % b = b;
        fval = fb;
        return
    elseif (fa > 0) == (fb > 0)
        error('MATLAB:fzero:ValuesAtEndPtsSameSign',...
            'The function values at the interval endpoints must differ in sign.')
    end
    
    % Starting guess scalar input
elseif (length(x) == 1)
    %%%%%%%%%%%%%%%%%%%%
    fx = fun(x,Y,dr);
    %%%%%%%%%%%%%%%%%%%%
    fcount = fcount + 1;  
    if fx == 0
        b = x;
        fval = fx;
        return
%     elseif ~isfinite(fx) || ~isreal(fx)
%         error('MATLAB:fzero:ValueAtInitGuessComplexOrNotFinite',...
%             'Function value at starting guess must be finite and real.');
    end
    
    if x ~= 0, 
        dx = x/50;
    else 
        dx = 1/50;
    end
    
    % Find change of sign.
    twosqrt = sqrt(2); 
    a = x; fa = fx; b = x; fb = fx;
    
    while (fa > 0) == (fb > 0)
        intervaliter = intervaliter + 1;
        dx = twosqrt*dx;
        a = x - dx;  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fa = fun(a,Y,dr);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fcount = fcount + 1;
        if (fa > 0) ~= (fb > 0) % check for different sign
            % Before we exit the while loop, print out the latest interval
            break
        end
        b = x + dx;  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fb = fun(b,Y,dr);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fcount = fcount + 1;        
    end % while
    savea = a; savefa = fa; saveb = b; savefb = fb;
else
    error('MATLAB:fzero:LengthArg2', 'Second argument must be of length 1 or 2.');
end % if (length(x) == 2

fc = fb;

% Main loop, exit from middle of the loop
while fb ~= 0 && a ~= b
    % Insure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
        c = a;  
        fc = fa;
        d = b - a;  
        e = d;
    end
    if abs(fc) < abs(fb)
        a = b;    
        b = c;    
        c = a;
        fa = fb;  
        fb = fc;  
        fc = fa;
    end
    
    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if (abs(m) <= toler) || (fb == 0.0) 
        break
    end
    
    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
        % Bisection
        d = m;  
        e = m;
    else
        % Interpolation
        s = fb/fa;
        if (a == c)
            % Linear interpolation
            p = 2.0*m*s;
            q = 1.0 - s;
        else
            % Inverse quadratic interpolation
            q = fa/fc;
            r = fb/fc;
            p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
            q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        end;
        if p > 0, 
            q = -q; 
        else
            p = -p; 
        end;
        % Is interpolated point acceptable
        if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))
            e = d;  
            d = p/q;
        else
            d = m;  
            e = m;
        end;
    end % Interpolation
    
    % Next point
    a = b;
    fa = fb;
    if abs(d) > toler, 
        b = b + d;
    elseif b > c, 
        b = b - toler;
    else
        b = b + toler;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fb = fun(b,Y,dr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fcount = fcount + 1;
    iter = iter + 1;
end % Main loop

fval = fb; % b is the best value



