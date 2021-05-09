% compute 余弦距离

function d = L2_distance_2(a,b)
% a,b: two matrices. each column is a data
% d:   distance matrix of a and b

a=a';b=b';%需要转置成n*d形式
num=size(a,1);
d=zeros(num,num);
parfor i=1:num
    for j=1:num
        d(i,j)=sum(a(i,:).*b(j,:))/((norm(a(i,:))*norm(b(j,:)))+eps);
    end
end
d=1-d;%
d = real(d);
%将距离规范化，可尝试多种规范化方法
    for  i = 1:num
        d(i,:) = (d(i,:) - min( d(i,:) ) ) / (max( d(i,:) )-min( d(i,:) )+eps) ;%eps防止分母为0
    end