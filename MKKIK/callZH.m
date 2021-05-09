function ZH = callZH(T,K,A0)

tau = size(A0,1);
nbkernel = size(K,3);
num = size(K,1);
TT = T*T';
ZH = zeros(nbkernel,1);
for p = 1 : nbkernel 
    tmp = 0;
    for i = 1 : num
      tmp = tmp + (1/num)*trace(K(A0(:,i),A0(:,i),p)*(eye(tau)-TT(A0(:,i),A0(:,i))));
    end
     ZH(p) = tmp;
end