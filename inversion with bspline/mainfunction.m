startpointx=[1,2,3,4,5,6];
startpointy=[1,2,3,4,5,6];
endpointx=[10,11,12,13];
endpointy=[10,11,12,13];
% 上面分别存储着起始点的x坐标和y坐标，这里只是个例子。

%存储数据个数
startnum=length(startpointx);
endnum=length(endpointx);

%如果想偷懒可以直接用这个部分的代码来确定范围，当然也可以自己来定义
xmin=floor(min(min(startpointx),min(endpointx)));
xmax=ceil(max(max(startpointx),max(endpointx)));

ymin=floor(min(min(startpointy),min(endpointy)));
ymax=ceil(max(max(startpointy),max(endpointy)));
square=1 ;%网格之间的距离

xnum=(xmax-xmin)/square;      %数据范围
ynum=(ymax-ymin)/square;
%G矩阵为对应的矩阵大小，G_raypath代表了路径密度的那个矩阵
G=zeros(startnum*endnum,xnum*ynum);
G_raypath=zeros(ynum,xnum);
for i=1:startnum
    for j=1:endnum
        xpoint=[startpointx(i),endpointx(j)];
        ypoint=[startpointy(i),endpointy(j)];
        [G1,G_raypath1]=getG(xpoint,ypoint,square,xmin,xmax,ymin.ymax);
             G((i-1)*endnum+j,:)=G1;
             G_raypath=G_raypath+G_raypath_1;
    end
end

%% 准备工作完成，下面进行计算
%由于基本的理论还需要测试完成，因此这部分只是一个初步的示例，并不代表最终的结果
%这部分的目的即把较小的敏感度矩阵的数值利用插值来进行替代，需要好好研究一下双线性插值这个算法
 G_sensitivity(1,i)=zeros(1,xum*ynum);
for i=1:xnum*ynum
    G_sensitivity(1,i)=sum(G(:,i));
end

G_sensitivity_matrix=reshape(G_sensitivity,ynum,xnum);


%% 在上面对G和M矩阵进行重建后，计算m的值




%% 多种情况下的误差分析