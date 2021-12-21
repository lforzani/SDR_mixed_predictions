function [A,B] = getABrr(C,d,invV)

maxiter = 1000;

[A,B] = initAB(C,d);

for iter = 1:maxiter
    B = updateB(A,invV,C);
    A = updateA(B,invV,C);
end


function [Aini,Bini] = initAB(C,d)
% esta funcion calcula A y B iniciales utilizando el estimador
% Dub-d

[U,D,V] = svd(C);

Ud = U(:,1:d);
Vd = V(:,1:d); 
Dd = D(1:d,1:d);

Aini = Ud*Dd;
Bini = Vd';


function B = updateB(A,invV,C)
% esta funcion acutaliza B utilizando el metodo iterativo
% de minimizacion cuadratica presentado en la pagina 47
% de la tesis de sabri

r = size(C,2);
d = size(A,2); 

kronn = kron(eye(r),A);
aux = inv(kronn'*invV*kronn)*kronn'*invV*C(:);

B = reshape(aux,[d,r]);


function A = updateA(B,invV,C)
% esta funcion acutaliza A utilizando el metodo iterativo
% de minimizacion cuadratica presentado en la pagina 47
% de la tesis de sabri
%Ojo que en la formula de la tesis de Sabri (47) la identidad del kron para
%el comupto de A esta mal. En vez de I_r es I_p2 (p2=m en nuestro caso).

m = size(C,1);
d = size(B,1); 
r = size(C,2);

kronn = kron(B, eye(m));
aux = inv(kronn*invV*kronn')*kronn*invV*C(:);

A = reshape(aux,[m,d]);




