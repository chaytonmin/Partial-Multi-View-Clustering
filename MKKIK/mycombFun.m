function cF = mycombFun(Y,gamma)

m = size(Y,3);
n = size(Y,1);
cF = zeros(n);
for p =1:m
    cF = cF + Y(:,:,p)*gamma(p);
end