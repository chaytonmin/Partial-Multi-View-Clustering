function h = fullLack(S,numpairedInst,singledNumView,num,v)
%补全缺失部分的S，即第二部分

h=zeros(num,num);

if v==2
  %取出S的存在的部分
    %缺失与共有，行
    s10=S(numpairedInst+1:numpairedInst+singledNumView,1:numpairedInst);
    s20=S(numpairedInst+singledNumView+1:num,1:numpairedInst);
   
    %共有与缺失，列
    s01=S(1:numpairedInst,numpairedInst+1:numpairedInst+singledNumView);
    s02=S(1:numpairedInst,numpairedInst+singledNumView+1:num);
   
    %计算需要补全缺失的S
    h12=s10*s02;
    h21=s20*s01;
    
    h(numpairedInst+1:numpairedInst+singledNumView,numpairedInst+singledNumView+1:num)=h12;
    h(numpairedInst+singledNumView+1:num,numpairedInst+1:numpairedInst+singledNumView)=h21;
end    

if v==3
    s10=S(numpairedInst+1:numpairedInst+singledNumView,1:numpairedInst);
    s20=S(numpairedInst+singledNumView+1:numpairedInst+singledNumView*2,1:numpairedInst);
    s30=S(numpairedInst+singledNumView*2+1:num,1:numpairedInst);
   
    %共有与缺失，列
    s01=S(1:numpairedInst,numpairedInst+1:numpairedInst+singledNumView);
    s02=S(1:numpairedInst,numpairedInst+singledNumView+1:numpairedInst+singledNumView*2);
    s03=S(1:numpairedInst,numpairedInst+singledNumView*2+1:num);
   
    %计算需要补全缺失的S
    h12=s10*s02;
    h13=s10*s03;
    h21=s20*s01;
    h23=s20*s03;
    h31=s30*s01;
    h32=s30*s02;
    
    h(numpairedInst+1:numpairedInst+singledNumView,numpairedInst+singledNumView+1:numpairedInst+singledNumView*2)=h12;
    h(numpairedInst+1:numpairedInst+singledNumView,numpairedInst+singledNumView*2+1:num)=h13;
    h(numpairedInst+singledNumView+1:numpairedInst+singledNumView*2,numpairedInst+1:numpairedInst+singledNumView)=h21;
    h(numpairedInst+singledNumView+1:numpairedInst+singledNumView*2,numpairedInst+singledNumView*2+1:num)=h23;
    h(numpairedInst+singledNumView*2+1:num,numpairedInst+1:numpairedInst+singledNumView)=h31;
    h(numpairedInst+singledNumView*2+1:num,numpairedInst+singledNumView+1:numpairedInst+singledNumView*2)=h32;
end

if v==6
    %取出S的存在的部分
    %缺失与共有，行
    s10=S(numpairedInst+1:numpairedInst+singledNumView,1:numpairedInst);
    s20=S(numpairedInst+singledNumView+1:numpairedInst+singledNumView*2,1:numpairedInst);
    s30=S(numpairedInst+singledNumView*2+1:numpairedInst+singledNumView*3,1:numpairedInst);
    s40=S(numpairedInst+singledNumView*3+1:numpairedInst+singledNumView*4,1:numpairedInst);
    s50=S(numpairedInst+singledNumView*4+1:numpairedInst+singledNumView*5,1:numpairedInst);
    s60=S(numpairedInst+singledNumView*5+1:num,1:numpairedInst);
    %共有与缺失，列
    s01=S(1:numpairedInst,numpairedInst+1:numpairedInst+singledNumView);
    s02=S(1:numpairedInst,numpairedInst+singledNumView+1:numpairedInst+singledNumView*2);
    s03=S(1:numpairedInst,numpairedInst+singledNumView*2+1:numpairedInst+singledNumView*3);
    s04=S(1:numpairedInst,numpairedInst+singledNumView*3+1:numpairedInst+singledNumView*4);
    s05=S(1:numpairedInst,numpairedInst+singledNumView*4+1:numpairedInst+singledNumView*5);
    s06=S(1:numpairedInst,numpairedInst+singledNumView*5+1:num);
    
    %计算需要补全缺失的S
    s12=s10*s02;s13=s10*s03;s14=s10*s04;s15=s10*s05;s16=s10*s06;
    s21=s20*s01;s23=s20*s03;s24=s20*s04;s25=s20*s05;s26=s20*s06;
    s31=s30*s01;s32=s30*s02;s34=s30*s04;s35=s30*s05;s36=s30*s06;
    s41=s40*s01;s42=s40*s02;s43=s40*s03;s45=s40*s05;s46=s40*s06;
    s51=s50*s01;s52=s50*s02;s53=s50*s03;s54=s50*s04;s56=s50*s06;
    s61=s60*s01;s62=s60*s02;s63=s60*s03;s64=s60*s04;s65=s60*s05;
    
    %补全缺失的S
    S(numpairedInst+1:numpairedInst+singledNumView,numpairedInst+singledNumView+1:numpairedInst+singledNumView*2)=s12;
    S(numpairedInst+1:numpairedInst+singledNumView,numpairedInst+singledNumView*2+1:numpairedInst+singledNumView*3)=s13;
    S(numpairedInst+1:numpairedInst+singledNumView,numpairedInst+singledNumView*3+1:numpairedInst+singledNumView*4)=s14;
    S(numpairedInst+1:numpairedInst+singledNumView,numpairedInst+singledNumView*4+1:numpairedInst+singledNumView*5)=s15;
    S(numpairedInst+1:numpairedInst+singledNumView,numpairedInst+singledNumView*5+1:num)=s16;
    
    
    S(numpairedInst+singledNumView+1:numpairedInst+singledNumView*2,numpairedInst+1:numpairedInst+singledNumView)=s21;
    S(numpairedInst+singledNumView+1:numpairedInst+singledNumView*2,numpairedInst+singledNumView*2+1:numpairedInst+singledNumView*3)=s23;
    S(numpairedInst+singledNumView+1:numpairedInst+singledNumView*2,numpairedInst+singledNumView*3+1:numpairedInst+singledNumView*4)=s24;
    S(numpairedInst+singledNumView+1:numpairedInst+singledNumView*2,numpairedInst+singledNumView*4+1:numpairedInst+singledNumView*5)=s25;
    S(numpairedInst+singledNumView+1:numpairedInst+singledNumView*2,numpairedInst+singledNumView*5+1:num)=s26;
    
    
    S(numpairedInst+singledNumView*2+1:numpairedInst+singledNumView*3,numpairedInst+1:numpairedInst+singledNumView)=s31;
    S(numpairedInst+singledNumView*2+1:numpairedInst+singledNumView*3,numpairedInst+singledNumView+1:numpairedInst+singledNumView*2)=s32;
    S(numpairedInst+singledNumView*2+1:numpairedInst+singledNumView*3,numpairedInst+singledNumView*3+1:numpairedInst+singledNumView*4)=s34;
    S(numpairedInst+singledNumView*2+1:numpairedInst+singledNumView*3,numpairedInst+singledNumView*4+1:numpairedInst+singledNumView*5)=s35;
    S(numpairedInst+singledNumView*2+1:numpairedInst+singledNumView*3,numpairedInst+singledNumView*5+1:num)=s36;
    
    S(numpairedInst+singledNumView*3+1:numpairedInst+singledNumView*4,numpairedInst+1:numpairedInst+singledNumView)=s41;
    S(numpairedInst+singledNumView*3+1:numpairedInst+singledNumView*4,numpairedInst+singledNumView+1:numpairedInst+singledNumView*2)=s42;
    S(numpairedInst+singledNumView*3+1:numpairedInst+singledNumView*4,numpairedInst+singledNumView*2+1:numpairedInst+singledNumView*3)=s43;
    S(numpairedInst+singledNumView*3+1:numpairedInst+singledNumView*4,numpairedInst+singledNumView*4+1:numpairedInst+singledNumView*5)=s45;
    S(numpairedInst+singledNumView*3+1:numpairedInst+singledNumView*4,numpairedInst+singledNumView*5+1:num)=s46;
    
    S(numpairedInst+singledNumView*4+1:numpairedInst+singledNumView*5,numpairedInst+1:numpairedInst+singledNumView)=s51;
    S(numpairedInst+singledNumView*4+1:numpairedInst+singledNumView*5,numpairedInst+singledNumView+1:numpairedInst+singledNumView*2)=s52;
    S(numpairedInst+singledNumView*4+1:numpairedInst+singledNumView*5,numpairedInst+singledNumView*2+1:numpairedInst+singledNumView*3)=s53;
    S(numpairedInst+singledNumView*4+1:numpairedInst+singledNumView*5,numpairedInst+singledNumView*3+1:numpairedInst+singledNumView*4)=s54;
    S(numpairedInst+singledNumView*4+1:numpairedInst+singledNumView*5,numpairedInst+singledNumView*5+1:num)=s56;
    
    S(numpairedInst+singledNumView*5+1:num,numpairedInst+1:numpairedInst+singledNumView)=s61;
    S(numpairedInst+singledNumView*5+1:num,numpairedInst+singledNumView+1:numpairedInst+singledNumView*2)=s62;
    S(numpairedInst+singledNumView*5+1:num,numpairedInst+singledNumView*2+1:numpairedInst+singledNumView*3)=s63;
    S(numpairedInst+singledNumView*5+1:num,numpairedInst+singledNumView*3+1:numpairedInst+singledNumView*4)=s64;
    S(numpairedInst+singledNumView*5+1:num,numpairedInst+singledNumView*4+1:numpairedInst+singledNumView*5)=s65;
end