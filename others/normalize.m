function X =normalize(X,v)
%%规范化方式
%%方法一将数据转化为以0为均值，既有正又有负数
%%方法二将数据转化到（0,1）之间
% %规范化方式一：z-score
num=size(X,1);

    for  j = 1:num
        X(j,:) = ( X(j,:) - mean( X(j,:) ) ) / (std( X(j,:) )+eps) ;%eps防止分母为0
    end

% %规范化方式二：max-min
% parfor i = 1 :v
%     for  j = 1:num
%         X{i}(j,:) = ( X{i}(j,:) - min( X{i}(j,:) ) ) / (max( X{i}(j,:) )-min( X{i}(j,:) )+eps) ;%eps防止分母为0
%     end
% end