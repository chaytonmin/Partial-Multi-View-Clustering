function X =normalize(X,v,num)
%%规范化方式
%%方法一将数据转化为以0为均值，既有正又有负数
%%方法二将数据转化到（0,1）之间
% %规范化方式一：z-score
for i = 1 :v
    for  j = 1:num
        X{i}(j,:) = ( X{i}(j,:) - mean( X{i}(j,:) ) ) / (std( X{i}(j,:) )+eps) ;%eps防止分母为0
    end
end
% %规范化方式二：max-min
% for i = 1 :v
%     for  j = 1:num
%         X{i}(j,:) = ( X{i}(j,:) - min( X{i}(j,:) ) ) / (max( X{i}(j,:) )-min( X{i}(j,:) )+eps) ;%eps防止分母为0
%     end
% end