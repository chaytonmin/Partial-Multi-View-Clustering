function [U,Yall,obj,acc_sym] = unified_subspace(L, Xc, X, Xa,gamma, beta,maxiter, Yc0, Y0,U0,gnd_all)
% Xc : complete views
%X  single
% Xa all 
Yc = Yc0;
U = U0;
Y = Y0;
k=size(Yc,2);
m=length(X);
Ik=eye(k,k);

for i=1:m
    d(i)=size(Xc{i},1);
    XX{i}=Xa{i}*Xa{i}';
    for j=1:m
        XLX{i,j}=Xa{i}*L{i,j}*Xa{j}';
    end
end 
maxiter=10;
for iter = 1:maxiter
    
    % Yc
    Apa=zeros(size(Yc));
    Ana=zeros(size(Yc));
    Gamma=zeros(size(Ik));
    for i=1:m
        A{i}=Xc{i}'*U{i};
        Ap{i}=(abs(A{i})+A{i})/2;
        Apa=Apa+Ap{i};
        An{i}=(abs(A{i})-A{i})/2;
        Ana=Ana+An{i}+Yc;
        Gammapart{i}=Yc'*A{i}-Ik;
        Gamma=Gamma+Gammapart{i};
    end
    Gammap=(abs(Gamma)+Gamma)/2;
    Gamman=(abs(Gamma)-Gamma)/2;
    Ycp=sqrt(Apa+Yc*Gamman)./(Ana+Yc*Gammap);
    Yc=Yc.*Ycp;
    
    %Y
   parfor i=1:m
        Y{i}=max(X{i}'*U{i},0);
        Ya{i}=[Yc;Y{i}];
    end
    
    %U
    for i=1:m
        Dii=2*sqrt(sum(U{i}.*U{i},2))+eps;
        D{i}=diag(1./Dii);
        %p1=inv(XX{i}+beta*D{i}+gamma*XLX{i,i});
        p1=(XX{i}+beta*D{i}+gamma*XLX{i,i});
        p2{i}=zeros(d(i),k);
        for j=1:m
            if i~=j
                p2{i}=p2{i}+XLX{i,j}*U{j};
            end
        end
        p3=Xa{i}*Ya{i}-gamma*p2{i};
        U{i}=p1\p3;
    end
     
     
     for i=1:m
     objp1(i)=norm(Xa{i}'*U{i}-Ya{i},'fro')^2;
     objp2(i)=sum(sqrt(sum(U{i}.*U{i},2)));
     for j=1:m
         obj3(i,j)=trace(U{i}'*Xa{i}*L{i,j}*Xa{j}'*U{j});
     end
     objp3=sum(obj3,2);
     end
     obj(iter)=sum(objp1+beta*objp2+gamma*objp3');
     
     Yall=Yc;
     for i=1:m
         Yall=[Yall;Y{i}];
     end
      
     [~,resultco]=max(Yall');
     
     [acc_sym(1,iter) acc_sym(2,iter)]=CalcMetrics(gnd_all, resultco'); % acc nmi
     acc_sym(3,iter)=RandIndex(gnd_all, resultco'); % ar 
     [acc_sym(4,iter) acc_sym(5,iter) acc_sym(6,iter)]= compute_f(gnd_all, resultco'); % f p r 

end

end
