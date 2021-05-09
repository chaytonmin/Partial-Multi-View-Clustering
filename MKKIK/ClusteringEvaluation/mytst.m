function [out,val] = mytst(a)
[inx1,inx2] = size(a);
avec = reshape(a,inx1*inx2,1);
indx = find(avec==max(avec));
indx1 = ceil(indx/inx1);
indx2 = mod(indx,inx1);
indx2(indx2==0)=inx1;
out=[indx2 indx1];
val = max(avec);