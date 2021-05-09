function [out,indx3] = mytst3D(a)

[len3] = size(a,3);
str = cell(len3,1);
val3 = zeros(len3,1);
for i3 =1:len3
    [str{i3}.out3, val3(i3)] = mytst(a(:,:,i3));
end
[val,indx] = max(val3);
indx3 = indx(end);
out = str{indx3}.out3;

