function [  ] = plotraypath1(startpointx,startpointy,endpointx,endpointy,figurenum)
%在给出了对应的起点和终点的数据组之后对数据进行绘制
%   前四个变量分别对应起始点的x,y坐标组以及终点的x,y坐标组，不对应的情况，由于这是使用模拟数据，因此暂时不考虑，最后一个是指的是figure的编号
figure(figurenum)
startnum=length(startpointx);
endnum=length(endpointx);
for i=1:startnum
    for j=1:endnum
        x=[startpointx(i),endpointx(j)];
        y=[startpointy(i),endpointy(j)];
        plot(x,y)
        hold on
    end
end
hold off
end

