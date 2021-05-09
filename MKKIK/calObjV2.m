function [obj,obj1,obj2]= calObjV2(T,K,gamma0)

nbkernel = size(K,3);
num = size(K,1);
ZH1 = zeros(nbkernel);
for p = 1 : nbkernel
    ZH1(p,p) = (trace(K(:,:,p))-trace(T'*K(:,:,p)*T));
end
obj = gamma0'*ZH1*gamma0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZH2 = zeros(nbkernel,1);
% for p =1:nbkernel
%     indexp = setdiff(1:nbkernel,p);
%     for q = 1:nbkernel-1
%         ZH2(p) = ZH2(p) + trace(K(:,:,p)*(-K(:,:,indexp(q))));
%         %% ZH2(p) = ZH2(p) + (1/(nbkernel*num))*(-trace(K(:,:,indexp(q))*K(:,:,p)'));
%     end
% %     for q = p:nbkernel
% %         %% ZH2(p) = ZH2(p) + (lambda/num)*trace(K(:,:,p)*(diag(sum(K(:,:,indexp(q))))-K(:,:,indexp(q))));
% %         ZH2(p,q) = (lambda/num)*trace(K(:,:,p)*K(:,:,q));
% %     end
% end
% % ZH2 = ZH2+ZH2'-diag(diag(ZH2));
% obj2 = gamma0'*ZH2;
% obj = obj1 + lambda*obj2;