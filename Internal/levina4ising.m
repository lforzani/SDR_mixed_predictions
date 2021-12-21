function gamma = levina4ising(H,z,lambda)
% INPUTS: 
% H: binarias, dimension nxq
% z: covariables, dimension pxq
% lambda:parametro de regularizacion para seleccion de variables. poner a
% cero inicialmente (valor por defecto).
%
% OUTPUTS:
% gamma: arreglo de coeficientes estimados. La dimension del arreglo
%        es (p+1) x q x q x length(lambda)
%
% =========================================================================
if nargin<3,
    lambda = 0;
end
p = size(z,2);
q = size(H,2);
L = length(lambda);
gamma_aux = SolPath(H,z,lambda,'joint');
gamma = zeros(p+1,q,q,L);
[gamma] = select_gamma(gamma_aux);
if length(lambda==1),
    gamma = squeeze(gamma);
end
    