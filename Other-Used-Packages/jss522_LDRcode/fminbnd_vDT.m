function [xf,fval] = fminbnd_vDT(fun,ax,bx,options,Y,dr)
% Descripcion minima:
% Inputs:
% - fun: handle a la función a optimizar. Dados Y y dr, 
% fun es una función solamente de x y se optimiza en el 
% intervalo [ax; bx]
% - options no lo usamos. Lo puse solamente por compatibilidad
% Outputs
% - xf: argumento que minimiza fun en el intervalo [ax,bx]
% - fval = fun(xf)



%%
% Definimos constantes para la terminar las iteraciones
tol = 1e-4;
maxiter = 500;
maxfun = 500; % este no se si tiene sentido...


funccount = 0;
% chequea si el limite superior para buscar la solución es mayor que el
% inferior. Sni no lo es, sale con error
if ax > bx,
    xf = []; fval = [];
    error('the upper bound must be larger than the lower bound');
end

% Computa un estimador inicial
seps = sqrt(eps);
c = 0.5*(3.0 - sqrt(5.0));
a = ax; 
b = bx;
v = a + c*(b-a);
w = v; 
xf = v;
d = 0.0; 
e = 0.0;
x= xf; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACA fun ES UN HANDLE A UNA FUNCION. 
% ESTO DEBERÍA TOMAR UNA FORMA DISTINTA DE ACUERDO AL MODELO USADO
fx = fun(x,Y,dr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
funccount = funccount + 1;

iter = 0;
fv = fx; 
fw = fx;
xm = 0.5*(a+b);
tol1 = seps*abs(xf) + tol/3.0;
tol2 = 2.0*tol1;


% Main loop
while ( abs(xf-xm) > (tol2 - 0.5*(b-a)) )
    gs = 1;
    % Is a parabolic fit possible?
    if abs(e) > tol1
        % Yes, so fit parabola
        gs = 0;
        r = (xf-w)*(fx-fv);
        q = (xf-v)*(fx-fw);
        p = (xf-v)*q-(xf-w)*r;
        q = 2.0*(q-r);
        if q > 0.0,  
            p = -p; 
        end
        q = abs(q);
        r = e;  
        e = d;

        % Is the parabola acceptable
        if ( (abs(p)<abs(0.5*q*r)) && (p>q*(a-xf)) && (p<q*(b-xf)) )
            % Yes, parabolic interpolation step
            d = p/q;
            x = xf+d;
            % f must not be evaluated too close to ax or bx
            if ((x-a) < tol2) || ((b-x) < tol2)
                si = sign(xm-xf) + ((xm-xf) == 0);
                d = tol1*si;
            end
        else
            % Not acceptable, must do a golden section step
            gs=1;
        end
    end
    if gs
        % A golden-section step is required
        if xf >= xm, 
            e = a-xf;    
        else
            e = b-xf;  
        end
        d = c*e;
    end

    % The function must not be evaluated too close to xf
    si = sign(d) + (d == 0);
    x = xf + si * max( abs(d), tol1 );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fu = fun(x,Y,dr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    funccount = funccount + 1;

    iter = iter + 1;

    % Update a, b, v, w, x, xm, tol1, tol2
    if fu <= fx
        if x >= xf, 
            a = xf; 
        else
            b = xf; 
        end
        v = w; 
        fv = fw;
        w = xf; 
        fw = fx;
        xf = x; 
        fx = fu;
    else % fu > fx
        if x < xf, 
            a = x; 
        else
            b = x; 
        end
        if ( (fu <= fw) || (w == xf) )
            v = w; fv = fw;
            w = x; fw = fu;
        elseif ( (fu <= fv) || (v == xf) || (v == w) )
            v = x; fv = fu;
        end
    end
    xm = 0.5*(a+b);
    tol1 = seps*abs(xf) + tol/3.0; 
    tol2 = 2.0*tol1;

    if funccount >= maxfun || iter >= maxiter
        fval = fx;
        return
    end
end % while
fval = fx;