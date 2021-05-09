function cF = mycombFun(Y,gamma)

m = size(Y,3); 
n1 = size(Y,1);
n2 = size(Y,2);
cF = zeros(n1,n2);
for p =1:m
    cF = cF + Y(:,:,p)*gamma(p);
end