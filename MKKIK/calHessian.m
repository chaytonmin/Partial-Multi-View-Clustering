function HE0 = calHessian(K,localFlag,A)

num = size(K,1);
numker = size(K,3);
HE0 = zeros(numker);
if localFlag==1
    for p =1:numker
        for q =p:numker
            tmp =0;
            for i =1:num
                tmp = tmp + (1/num)*trace(K(A(:,i),A(:,i),p)*K(A(:,i),A(:,i),q));
            end
            HE0(p,q) = tmp;
        end
    end
else
    for p =1:numker
        for q = p:numker
            HE0(p,q) = trace(K(:,:,p)'*K(:,:,q));
        end
    end
end
HE0 = (HE0+HE0')-diag(diag(HE0));