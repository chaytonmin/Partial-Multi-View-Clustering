% load('HWfold.mat');
% load('HWMean.mat');
% load('UCI_handwritten_digit.mat');
% s=X;

a=[4,3;5,6];b=[3,4;7,8];
num=2;
for i=1:num
    for j=1:num
        d(i,j)=sum(a(i,:).*b(j,:))/(norm(a(i,:))*norm(b(j,:)));
    end
end
d=1-d;%d:余弦相似度，1-d代表距离
d =normalize(d,1);%将距离规范化，可尝试多种规范化方法