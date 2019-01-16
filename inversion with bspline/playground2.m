%% 生成模型和模拟数据
% 2019.1.14需要补充完成ref情况下的解
clear
% x=10:10:100;
% x=x-5;
k1=1;
% 生成模拟数据方案1
% startnum=5;%这是以组来生成数据的，这里的是震源点每组startnum组，每组10个，接收点也类似
% endnum=3;
% y0interval=(1100-100)/(startnum-1);%这是之间的间隔
% y1interval=(1100-100)/(endnum-1);
% G=zeros(startnum*endnum*100,110);
% figure(1)
% for m=1:startnum
%     for n=1:endnum
%         for i=1:10
%             for j=1:10
% %                 ypoint=[y0interval*(m-1)+x(i),y1interval*(n-1)+x(j)];
%                 xpoint=[0,1000];
%                 plot(xpoint,ypoint)
%                 hold on
%                 [G1,G2]=getG(xpoint,ypoint,100,0,1000,0,1100);
%                 G(k1,:)=G1;
%                 k1=k1+1;
%             end
%         end
%     end
% end
% hold off

%生成模拟数据方案2
startnum=15;
endnum=5;
x0=0;
x1=1000;
y0=0;
y1=1000;
square=100;
interval1=(y1-y0)/(startnum);
interval2=(y1-y0)/(endnum);
xnum=(x1-x0)/square;
ynum=(y1-y0)/square;
G=zeros(startnum*endnum,(x1-x0)*(y1-y0)/(square*square));
G_raypath=zeros(ynum,xnum);
figure(1)
for m=1:startnum
    for n=1:endnum
        
%                 ypoint=[y0interval*(m-1)+x(i),y1interval*(n-1)+x(j)];
                ypoint=[(m-0.5)*interval1,(n-0.5)*interval2+0.01];
                xpoint=[x0,x1];
                plot(xpoint,ypoint)
                hold on
                [G1,G2]=getG(xpoint,ypoint,square,x0,x1,y0,y1); %完全在边缘点的时候会出现bug，暂时没有时间去修复这个bug了，自己注意一下
                G(k1,:)=G1;
                G_raypath=G_raypath+G2;
                k1=k1+1;
    end
end
title('raypath')
hold off
% for i=1:(y1-y0)/square
%     for j=1:(x1-x0)/square
%         m((y1-y0)/square*(i-1)+j,1)=3;
%     end
% end
AVE=mean(G_raypath);
G_raypathplot=ones(ynum+1,xnum+1);
G_raypathplot=G_raypathplot*10;
   figure(2)
  for i=1:ynum+1
     for j=1:xnum+1
         X1(i,j)=(j-1)*100;
         Y1(i,j)=(i-1)*100;
         if i<=ynum && j<=xnum
         G_raypathplot(i,j)=G_raypath(i,j);
         end
     end
  end
  %[XX,YY,ZZ]=griddata(X,Y,Z,linspace(min(X),max(X))',linspace(min(X),max(Y)),'v4');
%  [XX,YY]=meshgrid(X,Y);
%  ZZ=meshgrid(Z,Z);
 pcolor(X1,Y1,G_raypathplot)
%shading interp
grid on
 view(0,90)
 colorbar
 title('raypath density')
for i=1:(y1-y0)/square
    for j=1:(x1-x0)/square
        m((y1-y0)/square*(i-1)+j,1)=100+100*(-1)^(i+j);
    end
end
ref_c=mean(m);%参考速度，或者说是平均速度
d=G*m+random('norm',0,1000);
%% 插值的各种参数
%% 准备工作完成，下面进行计算
%由于基本的理论还需要测试完成，因此这部分只是一个初步的示例，并不代表最终的结果
%这部分的目的即把较小的敏感度矩阵的数值利用插值来进行替代，需要好好研究一下双线性插值这个算法


%重要参数！！
G_sensitivity_low=300;



G_sensitivity=zeros(1,xnum*ynum);
for i=1:xnum*ynum
    G_sensitivity(1,i)=sum(G(:,i));
end


%敏感度下限值

%将为0的矩阵元素设为-1
I=find(G_sensitivity<10e-15);
G_sensitivity(I)=-1;
%找到所有小点，设为0
[row,col]=find((G_sensitivity.*(G_sensitivity-G_sensitivity_low))<=10e-15);
G_sensitivity(row,col)=0;
G_sensitivity_matrix=reshape(G_sensitivity,xnum,ynum)';
clear G_sensitivity
G_sensitivity=G_sensitivity_matrix;
[row,col]=find((G_sensitivity_matrix.*(G_sensitivity_matrix-G_sensitivity_low))<=10e-15);
%% 进行检索并插值
small_num=length(row);
% point,pointa,,pointb,,pointc,,pointd,wa,wb,wc,wd     (x,y)
W_small=zeros(small_num,14);
%插值的距离之和上限
updistance=10;
nearupdistance=5; %临近插值的上限
for i=1:small_num
    W_small(i,1)=col(i);
    W_small(i,2)=row(i);% 记录求取的小点的坐标
    
    %a部分
    %datapointline 存储x_left,y_left,x_right,y_right,distance，bignum
    datapointline=zeros(ynum,6);
    %[datapointline(row(i),1:4),distance,bignum]=linesearch(G_sensitivity(row(i),:),row(i),row(i),col(i));
    [datapointline(row(i),1:4),datapointline(row(i),5),datapointline(row(i),6)]=linesearch(G_sensitivity(row(i),:),row(i),row(i),col(i));
    %b 部分
    datapointcol=zeros(xnum,6);
    [datapointcol(col(i),1:4),datapointcol(col(i),5),datapointcol(col(i),6)]=colsearch(G_sensitivity(:,col(i)),col(i),row(i),col(i));
    %c部分
    %c1 都得到两个点的情况
    if abs(datapointline(row(i),6)-2)<10e-15&&abs(datapointcol(col(i),6)-2)<10e-15
        if datapointline(row(i),5)<=datapointcol(col(i),5)
            [w1,w2]=lineinterp2(datapointline(row(i),1:4),col(i),row(i));
            W_small(i,3:6)=datapointline(row(i),1:4);
            W_small(i,11:12)=[w1,w2];
        else
            [w1,w2]=colinterp2(datapointcol(col(i),1:4),col(i),row(i));
            W_small(i,3:6)=datapointcol(col(i),1:4);
            W_small(i,11:12)=[w1,w2];
        end
        continue
    end
    %c2 只有一个2各点的情况
    if abs(datapointline(row(i),6)-2)<10e-15&&~(abs(datapointcol(col(i),6)-2)<10e-15)
        [w1,w2]=lineinterp2(datapointline(row(i),1:4),col(i),row(i));
        W_small(i,3:6)=datapointline(row(i),1:4);
        W_small(i,11:12)=[w1,w2];
        continue
    end
    if abs(~(datapointline(row(i),6)-2)<10e-15)&&(abs(datapointcol(col(i),6)-2)<10e-15)
        [w1,w2]=colinterp2(datapointcol(col(i),1:4),col(i),row(i));
        W_small(i,3:6)=datapointcol(col(i),1:4);
        W_small(i,11:12)=[w1,w2];
        continue
    end
    % c3 只有临近点的情况
    %c31 两个都有临界点的情况
    if abs(datapointline(row(i),6)+1)<10e-15&&abs(datapointcol(col(i),6)+1)<10e-15
        if datapointline(row(i),5)<=datapointcol(col(i),5)
            if datapointline(row(1),5)<nearupdistance
                W_small(i,3:4)=datapointline(row(i),1:2);
                W_small(i,11)=1;
                continue
            end
        else
            if datapointcol(col(i),5)<nearupdistance
                W_small(i,3:4)=datapointcol(col(i),1:2);
                W_small(i,11)=1;
                continue
            end
        end
    end
    %c31 只有一种临近点的情况
    if abs(datapointline(row(i),6)+1)<10e-15&&~(abs(datapointcol(col(i),6)+1)<10e-15)
        if datapointline(row(1),5)<nearupdistance
            W_small(i,3:4)=datapointline(row(i),1:2);
            W_small(i,11)=1;
            continue
        end
    end
    if ~abs((datapointline(row(i),6)+1)<10e-15)&&(abs(datapointcol(col(i),6)+1)<10e-15)
        if datapointcol(col(i),5)<nearupdistance
            W_small(i,3:4)=datapointcol(col(i),1:2);
            W_small(i,11)=1;
            continue
        end
        
    end
    % e部分
    %e1
    line_up=0;
    line_low=0;
    col_left=0;
    col_right=0;
    datapointlineup=inf(1,6);
    datapointlinelow=inf(1,6);
    datapointcolleft=inf(1,6);
    datapointcolright=inf(1,6);
    if row(i)-1>0
        for j=1:row(i)-1
            [datapointline(j,1:4),datapointline(j,5),datapointline(j,6)]=linesearch(G_sensitivity(j,:),j,row(i),col(i));
            if abs(datapointline(j,6)-2)<10e-15
                datapointline(j,5)=datapointline(j,5)./2;
            end
            if abs(datapointline(j,6)+1)<10e-15
                datapointline(j,5)=datapointline(j,5)*4;  %!!!!这是一种控制的方法，代表临近点插值并不完全受到欢迎，除非差的太远
            end
            if datapointline(j,6)>10e-15
                line_up=line_up+1;
            end
        end
        %筛选最小值点
        if line_up>10e-15
            [a,b]=min(datapointline(1:row(i)-1,5));
            datapointlineup=datapointline(b,:);
        end
    end
    %e2
    if row(i)<ynum %防止在最后一行
        for j=row(i)+1:ynum
            [datapointline(j,1:4),datapointline(j,5),datapointline(j,6)]=linesearch(G_sensitivity(j,:),j,row(i),col(i));
            if abs(datapointline(j,6)-2)<10e-15
                datapointline(j,5)=datapointline(j,5)./2;
            end
            if abs(datapointline(j,6)+1)<10e-15
                datapointline(j,5)=datapointline(j,5)*4;  %!!!!这是一种控制的方法，代表临近点插值并不完全受到欢迎，除非差的太远
            end
            if (datapointline(j,6))>10e-15
                line_low=line_low+1;
            end
        end
        %筛选最小值点
        if line_low>10e-15
            [a,b]=min(datapointline(row(i)+1:ynum,5));
            datapointlinelow=datapointline(b+row(i),:);
        end
    end
    
    %e3
    % f部分
    %f1
    if col(i)>0
        for j=1:col(i)-1
            [datapointcol(j,1:4),datapointcol(j,5),datapointcol(j,6)]=colsearch(G_sensitivity(:,j),j,row(i),col(i));
            if abs(datapointcol(j,6)-2)<10e-15
                datapointcol(j,5)=datapointcol(j,5)./2;
            end
            if abs(datapointcol(j,6)+1)<10e-15
                datapointcol(j,5)=datapointcol(j,5)*4;  %!!!!这是一种控制的方法，代表临近点插值并不完全受到欢迎，除非差的太远
            end
            if (datapointcol(j,6))>10e-15
                col_left=col_left+1;
            end
        end
        if col_left>10e-15
            [a,b]=min(datapointcol(1:col(i)-1,5));
            datapointcolleft=datapointcol(b,:);
        end
    end
    %筛选最小值点
    
    %f2  
    if col(i)<xnum
        for j=col(i)+1:xnum
            [datapointcol(j,1:4),datapointcol(j,5),datapointcol(j,6)]=colsearch(G_sensitivity(:,j),j,row(i),col(i));
            if abs(datapointcol(j,6)-2)<10e-15
                datapointcol(j,5)=datapointcol(j,5)./2;
            end
            if abs(datapointcol(j,6)+1)<10e-10
                datapointcol(j,5)
                datapointcol(j,5)=datapointcol(j,5)*4;  %!!!!这是一种控制的方法，代表临近点插值并不完全受到欢迎，除非差的太远
            end
            if (datapointcol(j,6))>10e-15
                
                col_right=col_right+1;
            end
            %选取最小值点
            if col_right>10e-15
                [a,b]=min(datapointcol(col(i)+1:xnum,5));
                datapointcolright=datapointcol(b+col(i),:);
            end
        end
    end
    %     %需要注意的是，如果想要加上距离的阈值判断的话，在这儿将line_up之类的参数设为0即可
    %     %暂时可以考虑不加因为我估计不会设置那么鬼畜的的路线分布
    %     %但如果有需要需要万分注意！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！加上数值
    %
    % h部分
    %h1
    if line_up*line_low*col_left*col_right>0
        if datapointlineup(1,5)+datapointlinelow(1,5)<=datapointcolleft(1,5)+datapointcolright(1,5)
            [w1,w2,w3,w4]=lineinterp4(datapointlineup(1,1:4),datapointlinelow(1,1:4),col(i),row(i));%注意这个函数里面需要判断一下哪个第二个元素为0，然后设置它的坐标和权重都是0
            W_small(i,3:10)=[datapointlineup(1,1:4),datapointlinelow(1,1:4)];
            W_small(i,11:14)=[w1,w2,w3,w4];
            continue
        else
            [w1,w2,w3,w4]=colinterp4(datapointcolleft(1,1:4),datapointcolright(1,1:4),col(i),row(i));%注意这个函数里面需要判断一下哪个第二个元素为0，然后设置它的坐标和权重都是0
            W_small(i,3:10)=[datapointcolleft(1,1:4),datapointcolright(1,1:4)];
            W_small(i,11:14)=[w1,w2,w3,w4];
            continue
        end
    end
    %h2
    if line_up*line_low>0 &&col_left*col_right<10e-16
        [w1,w2,w3,w4]=lineinterp4(datapointlineup(1,1:4),datapointlinelow(1,1:4),col(i),row(i));%注意这个函数里面需要判断一下哪个第二个元素为0，然后设置它的坐标和权重都是0
        W_small(i,3:10)=[datapointlineup(1,1:4),datapointlinelow(1,1:4)];
        W_small(i,11:14)=[w1,w2,w3,w4];
        continue
    end
    %     %h3
    if line_up*line_low<10e-16 &&col_left*col_right>0
        [w1,w2,w3,w4]=colinterp4(datapointcolleft(1,1:4),datapointcolright(1,1:4),col(i),row(i));%注意这个函数里面需要判断一下哪个第二个元素为0，然后设置它的坐标和权重都是0
        W_small(i,3:10)=[datapointcolleft(1,1:4),datapointcolright(1,1:4)];
        W_small(i,11:14)=[w1,w2,w3,w4];
        continue
    end
    %h4 检索插值得到的最近点(行或列插值得到)
    if line_up+line_low+col_left+col_right>0
        h4data=[datapointlineup;datapointlinelow;datapointcolleft;datapointcolright];
        [a,b]=min(h4data(:,5));
        h4updistance=6;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!上限值，非常重要
        if b>updistance
            continue
        end
        if b<=2
            [w1,w2]=lineinterp2(h4data(b,:),col(i),row(i));
            W_small(i,3:6)=h4data(b,1:4);
            W_small(i,11:12)=[w1,w2];
            continue
        else
            [w1,w2]=colinterp2(h4data(b,:),col(i),row(i));
            W_small(i,3:6)=h4data(b,1:4);
            W_small(i,11:12)=[w1,w2];
            continue
        end
    end
    %     %h5
    %     计算整个大点区域和这个点之间的差距，选择最小的那个，如果距离最小的且符合我们的要求，就是用最临近插值,因为在linesearch，colsearch里面已经加过限制了，因此在这儿就不加限制了里面已经加过
    for j=1:ynum
        if abs(datapointline(j,6)+1)<10e-16
            datapointline(j,5)=datapointline(j,5)/4;
        end
    end
    for j=1:xnum
        if abs(datapointcol(j,6)+1)<10e-16
            datapointcol(j,5)=datapointcol(j,5)/4;
        end
    end
    [a,b]=min(datapointline(:,5));
    h2data(1,:)=datapointline(b,:);
    [a,b]=min(datapointcol(:,5));
    h2data(2,:)=datapointcol(b,:);
    [a,b]=min(h2data(:,5));
    h2updistance=5;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!上限值，非常重要
    if a<=updistance
    W_small(i,3:4)=h2data(b,1:2);
    W_small(i,11)=1;
    end
 
 
end
 %% 利用权重矩阵W_small进行重组（注意这只是记录了信息，并不是新矩阵w）
 %更新：我已经完成了权重矩阵W的构成
% point,pointa,,pointb,,pointc,,pointd,wa,wb,wc,wd     (x,y)
G_new=G;
d_new=d;
for k=1:small_num
    if sum(W_small(k,11:14))<10e-10
        d_new=d_new-ref_c*G(:,W_small(k,1)+(W_small(k,2)-1)*xnum);
        continue
    end
    for t=1:4
        if (W_small(k,2+2*t))>0
    G_new(:,W_small(k,1+2*t)+(W_small(k,2+2*t)-1)*xnum)=G_new(:,W_small(k,1+2*t)+(W_small(k,2+2*t)-1)*xnum)+G(:,W_small(k,1)+(W_small(k,2)-1)*xnum)*W_small(k,10+t);
    G_new(:,W_small(k,1)+(W_small(k,2)-1)*xnum)=0;
        end
    end
end
[row_small,col_small]=find((G_sensitivity.*(G_sensitivity-G_sensitivity_low))<=10e-15);
[row_big,col_big]=find(((G_sensitivity-G_sensitivity_low))>10e-15);
nn=length(row_small);
W=eye(xnum*ynum,xnum*ynum);
for k=1:small_num
    if sum(W_small(k,11:14))<10e-10
        d_new=d_new-ref_c*G(:,W_small(k,1)+(W_small(k,2)-1)*xnum);
        W(W_small(k,1)+(W_small(k,2)-1)*xnum,W_small(k,1)+(W_small(k,2)-1)*xnum)=0;
        continue
    end
    for t=1:4
        if (W_small(k,2+2*t))>0
    W(W_small(k,1)+(W_small(k,2)-1)*xnum,W_small(k,1+2*t)+(W_small(k,2+2*t)-1)*xnum)=W_small(k,10+t);
    W(W_small(k,1)+(W_small(k,2)-1)*xnum,W_small(k,1)+(W_small(k,2)-1)*xnum)=0;
        end
    end
end
%% 储存大点小点和没有数据的点（零点）
%x,y,store_num,type（1代表大点，0代表小点，-1代表零点或者未能够进行插值的点）
m_store=zeros(ynum*xnum,1); %参考速度
num_big=0;
num_sma=0;
num_ref=0;
for i=1:ynum
    for j=1:xnum
        m_store(ynum*(i-1)+j,1)=j;
        m_store(ynum*(i-1)+j,2)=i;
        m_store(ynum*(i-1)+j,3)=ynum*(i-1)+j;
        if G_sensitivity(i,j)>10e-16
            m_store(ynum*(i-1)+j,4)=1;
            num_big=num_big+1;
        else
            if G_sensitivity(i,j)<-0.00005        %没有数据点
                m_store(ynum*(i-1)+j,4)=-1;
                num_ref=num_ref+1;
            end
        end
        if abs(G_sensitivity(i,j))<0.00005 %为0的情况
            if sum(W_small(k,11:14))>10e-10%有进行插值
                m_store(ynum*(i-1)+j,4)=0;
                num_sma=num_sma+1;
            else
                m_store(ynum*(i-1)+j,4)=-1;
                num_ref=num_ref+1;
            end
        end
    end
end
k1=1;
k2=1;
Gline=length(G(:,1));
%G_fresh=zeros(Gline,num_big+num_sma);
G_newfresh=zeros(Gline,num_big);
W_tmp=zeros(xnum*ynum,num_big);
W_fresh=zeros(xnum*ynum,num_big);
for n=1:xnum*ynum
    if sum(G(:,n))>1
     G_fresh(:,k1)=G(:,n);
     k1=k1+1;
    end
    if m_store(n,4)>0.1
     G_newfresh(:,k2)=G_new(:,n);
     k2=k2+1;
    end
end
% W清除计划
kcol=1;
kline=1;
for n=1:xnum*ynum    %清除列 %行没有清除的必要
    if m_store(n,4)>0.1
        W_fresh(:,kcol)=W(:,n);
        kcol=kcol+1;
    end
end
m_fresh=zeros(xnum*ynum,1);
for n=1:xnum*ynum    %清除列 %行没有清除的必要
    if m_store(n,4)<-0.1
       m_fresh(n,1)=ref_c;
    end
end
% for n=1:xnum*ynum    %清除行
%     if sum(W_tmp(n,:))>0.01
%         W_fresh(kline,:)=W_tmp(n,:);
%         kline=kline+1;
%     end
% end
%% 开始测试计算结果
% SVD


%  [U,S,V]=svd(G);
% % IU=find(abs(U)<10e-5);
% % U(IU)=0;
% % IS=find(abs(S)<10e-5);
% % S(IS)=0;
% % IV=find(abs(V)<10e-5);
% % V(IV)=0;
% n=85;
% for i=1:n
%     t(i)=abs(U(:,i)'*d/S(i,i));
% end
% figure(5)
% plot(t)
% Vp=V(:,1:n);
% Up=U(:,1:n);
% S_inverse=zeros(n,n);
% for i=1:n
%     S_inverse(i,i)=1/S(i,i);
% end
% Rm=diag(Vp*Vp');
% %m_USV=Vp*Vp'*m;
% m_USV=Vp*S_inverse*Up'*d;
% [U1,S1,V1]=svd(G_newfresh);
% % IU1=find(abs(U1)<10e-4);
% % U1(IU1)=0;
% % IS1=find(abs(S1)<10e-4);
% % S1(IS1)=0;
% % IV1=find(abs(V1)<10e-4);%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!注意abs
% % V1(IV1)=0;
% n=56;
% for i=1:n
%     t1(i)=abs(U1(:,i)'*d/S1(i,i));
% end
% figure(6)
% plot(t1)
% Vp1=V1(:,1:n);
% Up1=U1(:,1:n);
% S_inverse1=zeros(n,n);
% for i=1:n
%     S_inverse1(i,i)=1/S1(i,i);
% end
% Rm1=diag(Vp1*Vp1');
% 
% k=1;
% for i=1:xnum*ynum-num_ref
%     if m_store(i,4)==1
%         Rm_compare(k,1)=Rm(i,1);
%         m_compare(k,1)=m(i,1);
%         k=k+1;
%     end
% end
% 
% %m_U1VUS1=W_fresh*(Vp1*Vp1')*m_compare;
% m_U1VUS1=W_fresh*Vp1*S_inverse1*Up1'*d+m_fresh;
% figure(7)
% plot(Rm_compare,'r*')
% hold on
% plot(Rm1)
% figure(8)
% plot(m,'r-')
% hold on
% plot(m_USV,'g--')
% hold on
%  plot(m_U1VUS1,'b--')
%  hold off
% 
%  for i=1:ynum
%      for j=1:xnum
%          X((i-1)*xnum+j,1)=x0+(j-0.5)*square;
%          Y((i-1)*xnum+j,1)=y0+(i-0.5)*square;
%      end
%  end
%  %Z=m_U1VUS1;
%  Z=m_USV;
%  [XX,YY,ZZ]=griddata(X,Y,Z,linspace(min(X),max(X))',linspace(min(X),max(Y)),'v4');
% %  [XX,YY]=meshgrid(X,Y);
% %  ZZ=meshgrid(Z,Z);
%  figure(9)
%  pcolor(XX,YY,ZZ)
% shading interp
%  view(0,90)
%  colorbar
%  figure(10)
%  Z=m_U1VUS1;
%   [XX,YY,ZZ]=griddata(X,Y,Z,linspace(min(X),max(X))',linspace(min(X),max(Y)),'v4');
% %  [XX,YY]=meshgrid(X,Y);
% %  ZZ=meshgrid(Z,Z);
%  pcolor(XX,YY,ZZ)
% shading interp
%  view(0,90)
%  colorbar
%  figure(11)
%   Z=m;
%   [XX,YY,ZZ]=griddata(X,Y,Z,linspace(min(X),max(X))',linspace(min(X),max(Y)),'v4');
% %  [XX,YY]=meshgrid(X,Y);
% %  ZZ=meshgrid(Z,Z);
%  pcolor(XX,YY,ZZ)
% shading interp
%  view(0,90)
%  colorbar
 
%% test TR
%ONE 
L0=eye(num_big,num_big);
number=1000;
low=-5;
up=3;
interval=(up-low)/number;
t=low:interval:up;
%a_c_L=zeros(1000,1);
m_c_L0=zeros(number,3);
m_L0_value=zeros(num_big,number);
for i=1:number
    m_c_L0(i,1)=10^t(i);%α的数值
    m_L0_value(:,i)=getm(G_newfresh,d_new,L0,m_c_L0(i,1));
    m_c_L0(i,2)=norm(L0*m_L0_value(:,i));%||m||
    m_c_L0(i,3)=norm(G_newfresh* m_L0_value(:,i)-d_new); %||Gm-d||
end
%%
figure(20)
plot(log10(m_c_L0(:,3)),(m_c_L0(:,2)))
xlabel('||Gm-d||')
ylabel('||m||')
figure(21)
plot(m_c_L0(:,3),m_c_L0(:,2))
xlabel('||Gm-d||')
ylabel('||m||')
%xlim([0,0.005])
 m_L0_W=getm(G_newfresh,d_new,L0,10.185913880541170);

 m_L0_NEW=W_fresh*m_L0_W;
  figure(23)
 Z=m_L0_NEW;
  for i=1:ynum
     for j=1:xnum
         X((i-1)*xnum+j,1)=x0+(j-0.5)*square;
         Y((i-1)*xnum+j,1)=y0+(i-0.5)*square;
     end
 end
  [XX,YY,ZZ]=griddata(X,Y,Z,linspace(min(X),max(X))',linspace(min(X),max(Y)),'v4');
%  [XX,YY]=meshgrid(X,Y);
%  ZZ=meshgrid(Z,Z);
 pcolor(XX,YY,ZZ)
shading interp
 view(0,90)
 colorbar
 %%
% a=6.792036326171850e-04
% m=getm(G_newfresh,d_new,L0,6.792036326171850e-04);
% norm(L0*m)
% norm(G*m-d)
% figure(3)
% plot(x,m)
% xlabel('x')
% ylabel('m')
% title('a=6.792*10^-4')

%original
L0_or=eye(xnum*ynum,xnum*ynum);
number=1000;
low=-5;
up=3;
interval=(up-low)/number;
t=low:interval:up;
%a_c_L=zeros(1000,1);
m_c_L0_or=zeros(number,3);
m_L0_value_or=zeros(xnum*ynum,number);
for i=1:number
    m_c_L0_or(i,1)=10^t(i);%α的数值
    m_L0_value_or(:,i)=getm(G,d,L0_or,m_c_L0_or(i,1));
    m_c_L0_or(i,2)=norm(L0_or*m_L0_value_or(:,i));%||m||
    m_c_L0_or(i,3)=norm(G* m_L0_value_or(:,i)-d); %||Gm-d||
end
%%
m_L0_or=getm(G,d,L0_or,0.3311);%0.1076);

figure(24)
plot(log10(m_c_L0_or(:,3)),(m_c_L0_or(:,2)))
xlabel('||Gm-d||')
ylabel('||m||')
figure(25)
plot(m_c_L0_or(:,3),m_c_L0_or(:,2))
xlabel('||Gm-d||')
ylabel('||m||')
%xlim([0,0.005])

  figure(26)
 Z=m_L0_or;
  for i=1:ynum
     for j=1:xnum
         X((i-1)*xnum+j,1)=x0+(j-0.5)*square;
         Y((i-1)*xnum+j,1)=y0+(i-0.5)*square;
     end
 end
  [XX,YY,ZZ]=griddata(X,Y,Z,linspace(min(X),max(X))',linspace(min(X),max(Y)),'v4');
%  [XX,YY]=meshgrid(X,Y);
%  ZZ=meshgrid(Z,Z);
 pcolor(XX,YY,ZZ)
shading interp
 view(0,90)
 colorbar
   figure(27)
 Z=m;
  for i=1:ynum
     for j=1:xnum
         X((i-1)*xnum+j,1)=x0+(j-0.5)*square;
         Y((i-1)*xnum+j,1)=y0+(i-0.5)*square;
     end
 end
  [XX,YY,ZZ]=griddata(X,Y,Z,linspace(min(X),max(X))',linspace(min(X),max(Y)),'v4');
%  [XX,YY]=meshgrid(X,Y);
%  ZZ=meshgrid(Z,Z);
 pcolor(XX,YY,ZZ)
shading interp
 view(0,90)
 colorbar
 figure(28)
 Z=abs(m-m_L0_NEW);
  for i=1:ynum
     for j=1:xnum
         X((i-1)*xnum+j,1)=x0+(j-0.5)*square;
         Y((i-1)*xnum+j,1)=y0+(i-0.5)*square;
     end
 end
  [XX,YY,ZZ]=griddata(X,Y,Z,linspace(min(X),max(X))',linspace(min(X),max(Y)),'v4');
%  [XX,YY]=meshgrid(X,Y);
%  ZZ=meshgrid(Z,Z);
 pcolor(XX,YY,ZZ)
shading interp
caxis([0,150]);
 view(0,90)
 colorbar
 title('new way')
figure(29)
 Z=abs(m-m_L0_or);
  for i=1:ynum
     for j=1:xnum
         X((i-1)*xnum+j,1)=x0+(j-0.5)*square;
         Y((i-1)*xnum+j,1)=y0+(i-0.5)*square;
     end
 end
  [XX,YY,ZZ]=griddata(X,Y,Z,linspace(min(X),max(X))',linspace(min(X),max(Y)),'v4');
%  [XX,YY]=meshgrid(X,Y);
%  ZZ=meshgrid(Z,Z);
 pcolor(XX,YY,ZZ)
shading interp
caxis([0,150]);
 view(0,90)
 colorbar
 title('original way')
 
 %% 分辨率比较
  Rm_NEW=diag(Resolutionm(G_newfresh,L0,0.3311));
  Rm_or=diag(Resolutionm(G,L0_or,0.3311));
  k_Rm=1;
 for i=1:xnum*ynum
     if m_store(i,4)>0
         Rm_or_compare(k_Rm,1)=Rm_or(i,1);
         k_Rm=k_Rm+1;
     end
 end
 figure(30)
 plot(Rm_NEW-Rm_or_compare);
 title('Resolution of newway minus Resolution of original way')
   k_Rm=1;
   Rm_minus=zeros(ynum+1,xnum+1);
   Rm_minus=Rm_minus-0.1;
 for i=1:ynum
     for j=1:xnum
         if m_store(ynum*(i-1)+j,4)>0
             Rm_minus(i,j)=Rm_NEW(k_Rm,1)-Rm_or_compare(k_Rm,1);
             k_Rm=k_Rm+1;
         else
             Rm_minus(i,j)=-0.1;
         end
     end
 end
 figure(31)
 pcolor(X1,Y1,Rm_minus)
 colorbar