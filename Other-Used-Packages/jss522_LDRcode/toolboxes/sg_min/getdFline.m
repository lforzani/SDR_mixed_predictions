function df = getdFline(dF)
% DFLINE the derivative of the objective function, F, 
%	along the geodesic passing through Y in the direction H a distance t.
%
%	f= DFLINE(t,Y,H)
%	Y is expected to satisfy Y'*Y = I.
%	H is expected to satisfy H = tangent(Y,H0)
%	
% role	high level algorithm, basic window dressing for the fzero function
df = @evaluate_dFline;
    function dfval = evaluate_dFline(t,Y,H)
    	[Y,H] = move(Y,H,t);
        dfval = ip(Y,H,grad(dF,Y));
    end
end

