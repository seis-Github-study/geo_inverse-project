clear
x=10:10:100;
x=x-5;

k=1;
startnum=5;
endnum=3;
y0interval=(1100-100)/(startnum-1);
y1interval=(1100-100)/(endnum-1);
G=zeros(startnum*endnum*100,110);
figure(1)
for m=1:startnum
    for n=1:endnum
        for i=1:10
            for j=1:10
                ypoint=[y0interval*(m-1)+x(i),y1interval*(n-1)+x(j)];
                xpoint=[0,1000];
                plot(xpoint,ypoint)
                hold on
                [G1,G2]=getG(xpoint,ypoint,100,0,1000,0,1100);
                G(k,:)=G1;
                k=k+1;
            end
        end
    end
end
hold off
%%
[U,S,V]=svd(G);
Vp=V(:,1:86);
Rm=diag(Vp*Vp');
xnum=10;
ynum=11;
%% 准备工作完成，下面进行计算
%由于基本的理论还需要测试完成，因此这部分只是一个初步的示例，并不代表最终的结果
%这部分的目的即把较小的敏感度矩阵的数值利用插值来进行替代，需要好好研究一下双线性插值这个算法
G_sensitivity=zeros(1,xnum*ynum);
for i=1:xnum*ynum
    G_sensitivity(1,i)=sum(G(:,i));
end


%敏感度下限值
G_sensitivity_low=7000;
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
        [a,b]=min(hedata(:,5));
        h4updistance=5;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!上限值，非常重要
        if b>updistance
            break
        end
        if b<=2
            [w1,w2]=lineinterp2(h4data(b,:),col(i),row(i));
            W_small(i,3:6)=h4data(b,:);
            W_small(i,11:12)=[w1,w2];
            continue
        else
            [w1,w2]=colinterp2(h4data(b,:),col(i),row(i));
            W_small(i,3:6)=h4data(b,:);
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
    if a>updistance
        break
    end
    W_small(i,3:4)=h2data(b,1:2);
    W_small(i,11)=1;
    
    
end

%     %\h6
%     %没有任何举动，即代表最后权重矩阵里都是0，在重组G里面需要进行判断，如果这样的话！！！！！！！！！！！！！！！！！！！！！！！！
%     %则利用背景速度将其去除
%% 利用权重矩阵W_small进行重组（注意这只是记录了信息，并不是新矩阵w）
% point,pointa,,pointb,,pointc,,pointd,wa,wb,wc,wd     (x,y)
G_new=G;
%d_new=d;
for k=1:small_num
    if sum(W_small(k,11:14))<10e-10
        %d_new=d_new-ref_c*G(:,W_small(k,1)+(W_small(k,2)-1)*xnum);
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
        %d_new=d_new-ref_c*G(:,W_small(k,1)+(W_small(k,2)-1)*xnum);
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
% m_store=zeros(ynum*xnum,1);
% ref_c=3; %参考速度
% num_big=0;
% num_sma=0;
% num_ref=0;
% for i=1:ynum
%     for j=1:xnum
%         m_store(ynum*(i-1)+j,1)=j;
%         m_store(ynum*(i-1)+j,2)=i;
%         m_store(ynum*(i-1)+j,3)=ynum*(i-1)+j;
%         if G_sensitivity(i,j)>10e-16
%             m_store(ynum*(i-1)+j,4)=1;
%             num_big=num_big+1;
%         else
%             if G_sensitivity(i,j)<-0.00005        %没有数据点
%                 m_store(ynum*(i-1)+j,4)=-1;
%                 num_ref=num_ref+1;
%             end
%         end
%         if abs(G_sensitivity(i,j))<0.00005 %为0的情况
%             if sum(W_small(k,11:14))>10e-10%有进行插值
%                 m_store(ynum*(i-1)+j,4)=0;
%                 num_sma=num_sma+1;
%             else
%                 m_store(ynum*(i-1)+j,4)=-1;
%                 num_ref=num_ref+1;
%             end
%         end
%     end
% end
% k1=1;
% k2=1;
% G_fresh=zeros()
% for n=1:xnum*ynum
%     if sum(G(:,n))>10-16
%      G_fresh(:,k1)=G(:,n);
%      k1=k1+1;
%     end
% end
%     [U1,S1,V1]=svd(G_new);
%     for i=1:84
%         S1_inv(i,i)=1/S1(i,i);
%     end
%     Vp1=V1(:,1:84);
%     Up1=U1(:,1:84);
%     Rm1=diag(Vp1*S1_inv*Up1'*G);
