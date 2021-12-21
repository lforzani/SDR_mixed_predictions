function [resultado1, resultado2]=diferencia(a,b)

aux1=a*inv(a'*a)*a';
aux2=b*inv(b'*b)*b';

resultado1=norm(aux1-aux2,'fro');
           resultado2=norm(aux1-aux2,2);
