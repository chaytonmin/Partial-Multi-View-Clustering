function [result]= MLAN_Clustering(X,c,groundtruth,options)
% X:                cell array, 1 by view_num, each array is num by d_v
% c:                number of clusters
% v:                number of views
% k:                number of adaptive neighbours
%lambda             参数
%alpha              S的参数
%beta：             缺失项的参数
%groundtruth：      groundtruth of the data, num by 1
%pairPortion        完整比例

I=options.I;
J=options.J;
lambda=options.lambda;
k=options.k;
alpha=options.alpha;
pairPortion=options.pairPortion;
fai=options.fai;
VIR=options.VIR;
v= size(X,2);
num = size(X{1},1);
NITER = 30;%迭代的最高次数
pairedNum=floor(num*pairPortion);%共有的个数
singledNum=ceil((num-pairedNum)*VIR);%单有的个数
% c1=repmat((pairedNum+singledNum)/num,num,1);
c1=ones(num,1);%初始时缺失之和等于1
% c2=repmat((1-(pairedNum+singledNum)/num),num,1);
c2=ones(num,1);%初始时缺失之和等于1
%% =====================  Initialization =====================
% for i = 1 :v
%     for  j = 1:num
%         X{i}(j,:) = ( X{i}(j,:) - mean( X{i}(j,:) ) ) / (std( X{i}(j,:) )+eps) ;%eps防止分母为0
%     end
% end
%initialize weighted_distX
SUM = zeros(num);
%   d= I{3}.* L2_distance_1( X{3}',X{3}' );   %用于测试该值
%    for  j = 1:size(d,1)
%                  d(j,:) = (d(j,:) - min( d(j,:) ) ) / (max( d(j,:) )-min( d(j,:) )+eps) ;%eps防止分母为0
%    end
for i = 1:v
    %%此处需要根据view数据类型，手动选择计算距离方法，L2_distance_1是欧式距离，L2_distance_1是余弦距离，均已规范化
%       if i==1
          d= I{i}.* L2_distance_1( X{i}',X{i}' );%图片使用欧式距离L1，文字使用余弦距离L2，例如BDGP数据，数据2      
%      else
%           d= I{i}.* L2_distance_2( X{i}',X{i}' );
%       end
      for  j = 1:size(d,1)
                 d(j,:) = (d(j,:) - min( d(j,:) ) ) / (max( d(j,:) )-min( d(j,:) )+eps) ;%eps防止分母为0
      end
%          for  j = 1:size(d,1)
%                  d(j,:) = (d(j,:) - mean( d(j,:) ) ) / (std( d(j,:)) +eps) ;%eps防止分母为0
%          end
      distX_initial(:,:,i)=d;
    SUM = SUM + distX_initial(:,:,i);
end

% SUM=fullSUM(SUM,pairedNum,singledNum,num,v);%初始化补全距离
distX=(1/v)*J.*SUM;
[distXs, idx] = sort(distX,2);

%initialize S
S = zeros(num);
rr = zeros(num,1);
for i = 1:pairedNum
    di = distXs(i,2:k+2);
     rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx(i,2:k+2);
    S(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps); %initialize S，利用第k+1个的距离与该距离之差比上k个k+1的距离与前k个距离之和，eps是防止为0
end;
% for i = pairedNum+1:num
%         b=num-singledNum-pairedNum+2; 
%         di = distXs(i,b:k+b);
%         rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
%         id = idx(i,b:k+b);
%         S(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
% end
for i = pairedNum+1:pairedNum+singledNum
        b=num-singledNum-pairedNum+2; 
        di = distXs(i,b:k+b);
        rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
        id = idx(i,b:k+b);
        S(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end
for i = pairedNum+singledNum+1:num
        b=singledNum+2; 
        di = distXs(i,b:k+b);
        rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
        id = idx(i,b:k+b);
        S(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end
beta = mean(rr);%初始化beta

% initialize F
S = (S+S')/2;                                                        
D = diag(sum(S));
L = D - S;
[F, temp, evs]=eig1(L, c, 0);

if sum(evs(1:c+1)) < 0.00000000001
    error('The original graph has more than %d connected component', c);
end;
p=1;

%% =====================  updating =====================
for iter = 1:NITER
    % update weighted_distX
    SUM = zeros(num,num);
    for i = 1 : v
        if iter ==1  
            distX_updated = distX_initial;
        end
             Wv(i) = (0.5*p)/(sum(sum(I{i}.* distX_updated(:,:,i).*S))+eps)^((2-p)/2); %加入系数
        distX_updated(:,:,i) =Wv(i)*distX_updated(:,:,i) ;
        SUM = SUM + distX_updated(:,:,i);
    end
    distX = J.*SUM;
    [distXs, idx] = sort(distX,2);
    
    %update S
    %update S属于视图信息存在,即第一部分
    distf = L2_distance_1(F',F');
    S = zeros(num);
    
      for i = 1:pairedNum
        a=2;     
        idxa0 = idx(i,a:a+k-1);
        dfi = distf(i,idxa0);
        dxi = distX(i,idxa0);
        ad = -(dxi+lambda*dfi)/(2*beta);
        S(i,idxa0) = EProjSimplex_new(ad);   
      end
%       for i = pairedNum+1:num
%         b=num-singledNum-pairedNum+2; 
%         idxa0 = idx(i,b:b+k-1);
%         dfi = distf(i,idxa0);
%         dxi = distX(i,idxa0);
%         ad = -(dxi+lambda*dfi)/(2*beta);
%         S(i,idxa0) = EProjSimplex_new(ad,c1(i));
%         c22(i)=1-sum(S(i,idxa0));
%        end
       for i = pairedNum+1:pairedNum+singledNum
        b=num-singledNum-pairedNum+2; 
        idxa0 = idx(i,b:b+k-1);
        dfi = distf(i,idxa0);
        dxi = distX(i,idxa0);
        ad = -(dxi+lambda*dfi)/(2*beta);
        S(i,idxa0) = EProjSimplex_new(ad,c1(i));
        c22(i)=1-sum(S(i,idxa0));
       end
       for i = pairedNum+singledNum+1:num
        b=singledNum+2; 
        idxa0 = idx(i,b:b+k-1);
        dfi = distf(i,idxa0);
        dxi = distX(i,idxa0);
        ad = -(dxi+lambda*dfi)/(2*beta);
        S(i,idxa0) = EProjSimplex_new(ad,c1(i));
        c22(i)=1-sum(S(i,idxa0));
       end
%        
% %      update S 属于视图信息缺失，填充缺失的相似度
         h = fullLack(S,pairedNum,singledNum,num,v);%利用得到的相似度计算缺失的相似度

         for i=pairedNum+1:num
             sumS=0;
             for j=1:num
                 if h(i,j)~=0
                     S(i,j)=(fai-((lambda*distf(i,j)-2*alpha*h(i,j)))/(2*(alpha+beta)));
                     S(i,j);
                     sumS=sumS+S(i,j);
                 end                 
             end
             c1(i)=1-sumS;
         end
%         save('D:\文献代码\S','S');
    %update F
    S = (S+S')/2;                                                      
    D = diag(sum(S));
    L = D-S;
    F_old = F;
    [F, temp, ev]=eig1(L, c, 0);
    evs(:,iter+1) = ev;
    
    %画迭代图
%     A=distX.*S;
    %obj(iter)=trace(F'*L*F);
%     obj(iter)=sum(A(:))+beta*norm(S,'fro')+2*lambda*trace(F'*L*F);
    
    %update lambda  代码中方法
%     thre = 1*10^-10;
%     fn1 = sum(ev(1:c));                                                
%     fn2 = sum(ev(1:c+1));
%     if fn1 > thre
%         lambda = 2*lambda;   %初始乘2
%     elseif fn2 < thre
%         lambda = lambda/2;  F = F_old;
%     else
%         break;
%     end;
%  lambda
%  iter
     %update lambda  实际论文讲的方法
    [clusternum, y]=graphconncomp(sparse(S)); 
     if clusternum > c
        lambda = lambda/2;F = F_old;
    elseif clusternum < c
        lambda = 2*lambda; 
     else
        break;
    end;
% % lambda   
% % iter
   %sprintf('iter = %d',iter);
%    [clusternum, y]=graphconncomp(sparse(S)); y = y';
% 
%     result = ClusteringMeasure(groundtruth, y);
%     nmi=result(2)
%     obj2(iter)=nmi;
end
%% =====================  result =====================
% plot(obj);
[clusternum, y]=graphconncomp(sparse(S)); y = y';%一种聚类算法，此算法不用k-means
if clusternum ~= c
    sprintf('Can not find the correct cluster number: %d', c)
end;
result= ClusteringMeasure(groundtruth, y);