G1=[0,0,0,  0,0,0,       0,0,0       0,0,0   56,0.085,12,   0,5,54,     1,5,1,     4,5,6];
G2=[0,0,0,  0,0,0,       0,0,0       0,0,0   56,0.085,12,   0,0,54,     1,5,1,     4,5,6];
G3=[0,0,0,  0,0,0,       0,0,0       0,0,0   56,0.085,12,   0,0,54,     1,5,1,     4,5,6];
G4=[0,0,0,  0,0,0,       0,0,0       0,0,0   56,0.085,12,   0,0,54,     1,5,1,     4,5,6];
G5=[0,0,0,  0,0,0,       0,0,0       0,0,0   56,0.085,12,   0,0,54,     1,5,1,     4,5,6];
G=[G1;G2;G3;G4;G5];
xnum=3;
ynum=8;
%% 准备工作完成，下面进行计算
%由于基本的理论还需要测试完成，因此这部分只是一个初步的示例，并不代表最终的结果
%这部分的目的即把较小的敏感度矩阵的数值利用插值来进行替代，需要好好研究一下双线性插值这个算法
G_sensitivity=zeros(1,xnum*ynum);
for i=1:xnum*ynum
    G_sensitivity(1,i)=sum(G(:,i));
end


%敏感度下限值
G_sensitivity_low=2;
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
for i=1:1
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
%     %c部分
%     %c1 都得到两个点的情况
%     if abs(datapointline(row(i),6)-2)<10e-15&&abs(datapointcol(col(i),6)-2)<10e-15
%         if datapointline(row(i),5)<=datapointcol(col(i),5)
%             [w1.w2]=lineinterp2(datapointline(row(i),1:4),col(i),row(i));
%             W_small(i,3:6)=datapointline(row(i),1:4);
%             W_small(i,11:12)=[w1,w2];
%         else
%             [w1.w2]=colinterp2(datapointcol(col(i),1:4),col(i),row(i));
%             W_small(i,3:6)=datapointcol(col(i),1:4);
%             W_small(i,11:12)=[w1,w2];
%         end
%         continue
%     end
%     %c2 只有一个2各点的情况
%     if abs(datapointline(row(i),6)-2)<10e-15&&~(abs(datapointcol(col(i),6)-2)<10e-15)
%         [w1.w2]=lineinterp2(datapointline(row(i),1:4),col(i),row(i));
%         W_small(i,3:6)=datapointline(row(i),1:4);
%         W_small(i,11:12)=[w1,w2];
%         continue
%     end
%     if abs(~(datapointline(row(i),6)-2)<10e-15)&&(abs(datapointcol(col(i),6)-2)<10e-15)
%         [w1.w2]=colinterp2(datapointcol(col(i),1:4),col(i),row(i));
%         W_small(i,3:6)=datapointcol(col(i),1:4);
%         W_small(i,11:12)=[w1,w2];
%         continue
%     end
%     % c3 只有临近点的情况
%     %c31 两个都有临界点的情况
%     if abs(datapointline(row(i),6)+1)<10e-15&&abs(datapointcol(col(i),6)+1)<10e-15
%         if datapointline(row(i),5)<=datapointcol(col(i),5)
%             W_small(i,3:4)=datapointline(row(i),1:2);
%             W_small(i,11)=1;
%         else
%             W_small(i,3:4)=datapointcol(col(i),1:2);
%             W_small(i,11)=1;
%         end
%         continue
%     end
%     %c31 只有一种临近点的情况
%     if abs(datapointline(row(i),6)+1)<10e-15&&~(abs(datapointcol(col(i),6)+1)<10e-15)
%         W_small(i,3:4)=datapointline(row(i),1:2);
%         W_small(i,11)=1;
%         continue
%     end
%     if abs(~(datapointline(row(i),6)+1)<10e-15)&&(abs(datapointcol(col(i),6)+1)<10e-15)
%         W_small(i,3:4)=datapointcol(col(i),1:2);
%         W_small(i,11)=1;
%         continue
%     end
%     % e部分
%     %e1
%     line_up=0;
%     line_low=0;
%     col_left=0;
%     col_right=0;
%     datapointlineup=inf(1,6);
%     datapointlinelow=inf(1,6);
%     datapointcolleft=inf(1,6);
%     datapointcolright=inf(1,6);
%     for j=1:row(i)-1
%         [datapointline(j,1:4),datapointline(j,5),datapointline(j,6)]=linesearch(G_sensitivity(j,:),j,row(i),col(i));
%         if abs(datapointline(j,6)-2)<10e-15
%             datapointline(j,6)=datapointline(j,5)./2;
%         end
%         if datapointline(j,6)>10e-15
%             line_up=line_up+1;
%         end
%     end
%     %筛选最小值点
%     if line_up>10e-15
%         [a,b]=min(datapointline(1:row(i)-1,5));
%         datapointlineup=datapointline(b,:);
%     end
%     %e2
%     for j=row(i)+1:ynum
%         [datapointline(j,1:4),datapointline(j,5),datapointline(j,6)]=linesearch(G_sensitivity(j,:),j,row(i),col(i));
%         if abs(datapointline(j,6)-2)<10e-15
%             datapointline(j,5)=datapointline(j,5)./2;
%         end
%         if (datapointline(j,6))>10e-15
%             line_low=line_low+1;
%         end
%     end
%     %筛选最小值点
%     if line_low>10e-15
%         [a,b]=min(datapointline(row(i)+1:ynum,5));
%         datapointlinelow=datapointline(b,:);
%     end
%     %e3
%     % f部分
%     %f1
%     for j=1:col(i)-1
%         [datapointcol(j,1:4),datapointcol(j,5),datapointcol(j,6)]=colsearch(G_sensitivity(col(i),:),col(i),row(i),col(i));
%     end
%     if abs(datapointcol(j,6)-2)<10e-15
%         datapointcol(j,6)=datapointcol(j,5)./2;
%     end
%     if (datapointcol(j,6))>10e-15
%         col_left=col_left+1;
%     end
%     %筛选最小值点
%     if col_left>10e-15
%         [a,b]=min(datapointlcol(1:col(i)-1,5));
%         datapointcolleft=datapointcol(b,:);
%     end
%     %f2
%     for j=col(i)+1:xnum
%         [datapointcol(j,1:4),datapointcol(j,5),datapointcol(j,6)]=colsearch(G_sensitivity(col(i),:),col(i),row(i),col(i));
%     end
%     if abs(datapointcol(j,6)-2)<10e-15
%         datapointcol(j,5)=datapointcol(j,5)./2;
%     end
%     if (datapointcol(j,6))>10e-15
%         col_right=col_right+1;
%     end
%     %选取最小值点
%     if col_right>10e-15
%         [a,b]=min(datapointlcol(col(i)+1:xnum,5));
%         datapointcolleft=datapointcol(b,:);
%     end
%     
%     %需要注意的是，如果想要加上距离的阈值判断的话，在这儿将line_up之类的参数设为0即可
%     %暂时可以考虑不加因为我估计不会设置那么鬼畜的的路线分布
%     %但如果有需要需要万分注意！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！加上数值
%     
%     % h部分
%     %h1
%     if line_up*line_low*col_left*col_right>0
%         if datapointlineup(1,5)+datapointlinelow(1,5)<=datapointcolleft(1,5)+datapointright(1,5)
%             [w1,w2,w3,w4]=lineinterp4(datapointlineup(1,1:4),datapointlinelow(1,1:4),col(i),row(i));%注意这个函数里面需要判断一下哪个第二个元素为0，然后设置它的坐标和权重都是0
%             W_small(i,3:10)=[datapointlineup(1,1:4),datapointlinelow(1,1:4)];
%             W_small(i,11:14)=[w1,w2,w3,w4];
%             continue
%         else
%             [w1,w2,w3,w4]=colinterp4(datapointcolleft(1,1:4),datapointcolright(1,1:4),col(i),row(i));%注意这个函数里面需要判断一下哪个第二个元素为0，然后设置它的坐标和权重都是0
%             W_small(i,3:10)=[datapointcolleft(1,1:4),datapointcolright(1,1:4)];
%             W_small(i,11:14)=[w1,w2,w3,w4];
%             continue
%         end
%     end
%     %h2
%     if line_up*line_low>0 &&col_left*col_right<10e-16
%         [w1,w2,w3,w4]=lineinterp4(datapointlineup(1,1:4),datapointlinelow(1,1:4),col(i),row(i));%注意这个函数里面需要判断一下哪个第二个元素为0，然后设置它的坐标和权重都是0
%         W_small(i,3:10)=[datapointlineup(1,1:4),datapointlinelow(1,1:4)];
%         W_small(i,11:14)=[w1,w2,w3,w4];
%         continue
%     end
%     %h3
%     if line_up*line_low<10e-16 &&col_left*col_right>0
%         [w1,w2,w3,w4]=colinterp4(datapointcolleft(1,1:4),datapointcolright(1,1:4),col(i),row(i));%注意这个函数里面需要判断一下哪个第二个元素为0，然后设置它的坐标和权重都是0
%         W_small(i,3:10)=[datapointcolleft(1,1:4),datapointcolright(1,1:4)];
%         W_small(i,11:14)=[w1,w2,w3,w4];
%         continue
%     end
%     %h4
%     if line_up+line_low+col_left+col_right>0
%         h4data=[datapointlineup;datapointlinelow;datapointcolleft;datapointcolright];
%         [a,b]=min(hedata(:,5));
%         h4updistance=5;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!上限值，非常重要
%         if b>updistance
%             break
%         end
%         if b<=2
%             [w1,w2]=lineinterp2(h4data(b,:),col(i),row(i));
%             W_small(i,3:6)=h4data(b,:);
%             W_small(i,11:12)=[w1,w2];
%             continue
%         else
%             [w1,w2]=colinterp2(h4data(b,:),col(i),row(i));
%             W_small(i,3:6)=h4data(b,:);
%             W_small(i,11:12)=[w1,w2];
%             continue
%         end
%     end
%     %h5 计算整个大点区域和这个点之间的差距，选择最小的那个，如果距离最小的且符合我们的要求，就是用最临近插值
%     %\h6
%     %没有任何举动，即代表最后权重矩阵里都是0，在重组G里面需要进行判断，如果这样的话！！！！！！！！！！！！！！！！！！！！！！！！
%     %则利用背景速度将其去除
%     
 end