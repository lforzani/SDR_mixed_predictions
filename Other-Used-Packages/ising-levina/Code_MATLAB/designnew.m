
function[y_new x_new] = designnew(y,x,option)
n = size(y,1);
q = size(y,2);
p = size(x,2);

y_new = reshape(y,n*q,1);

K = zeros(n, (p+1)*q);
for i = 1:n
K(i,:) = kron(y(i,:),[1 x(i,:)]);
end

if strcmp(option,'separate')
    x_new = zeros(n*q, p+(p+1)*(q-1));
    x_new(:,1:p) = repmat(x,q,1);
    for j = 1:q
        tmp = K;
        tmp(:,(j-1)*(p+1)+1:j*(p+1))=[];
        x_new(n*(j-1)+1:n*j, p+1:end) = tmp;
    end
    
elseif strcmp(option, 'joint')
    x_new = zeros(n*q,q*(q+1)*(p+1)/2);
    z = zeros(n,p+1,q); 
  	for j = 1:q
  		z(:,:,j) = K(:,((j-1)*(p+1)+1):(j*(p+1)));
    end
  	x_new(:,1:(p+1)*q) = kron(eye(q),[ones(n,1),x]);
    for j = 1:q-1
        for k = j+1:q
            r = j*(2*q-j-1)/2+k;
            x_new(((j-1)*n+1):(j*n), ((r-1)*(p+1)+1):(r*(p+1))) = z(:,:,k);
     		x_new(((k-1)*n+1):(k*n), ((r-1)*(p+1)+1):(r*(p+1))) = z(:,:,j);
        end
    end
end

end
    
                
    
    
    

















