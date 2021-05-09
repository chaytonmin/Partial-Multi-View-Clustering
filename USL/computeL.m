function [L]=computeL(Xa, gnd_all,options) 

 views=length(Xa);
 W=cell(views,views);
 L=cell(views,views);

for i=1:views
    for j=1:views
       N1=length(gnd_all{1,i});
       N2=length(gnd_all{1, j});
        if i==j
    W{i,j}=constructW_cai(Xa{i},options);
        else
    W{i,j}= inter_view(gnd_all{i},gnd_all{j},N1,N2); 
        end
      DCol = full(sum(W{i,j},2));
    D = spdiags(DCol,0,N1,N2);
    L{i,j} = D - W{i,j};  
    end
end


function  Wtr=inter_view(gt1,gt2,n1,n2)
Wtr=zeros(n1,n2);
for i=1:n1
    for j=1:n2
       if gt1(i)==gt2(j)
           Wtr(i,j)=1;
       end
    end
end
