function [x,y,ix] = delp(b,st)

global p;
global M;
global N;

t = 0;
s = 0;
y = st;


for i=1:p
    
   if st(i) == 0  
       continue
   else
       t = t +1;
       if norm(b(t,:)) <= 1e-6
           s = s+1;
           y(i)=0;
           ix(s)=t;
       else
           continue
       end
   end
end

if sum(y) < sum(st)
b(ix,:)=[];
M(:,ix)=[];
M(ix,:)=[];
N(:,ix)=[];
N(ix,:)=[];
x = b; 
else
    ix=0;
    x=b;
end
