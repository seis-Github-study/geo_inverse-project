function [  ] = plotraypath1(startpointx,startpointy,endpointx,endpointy,figurenum)
%�ڸ����˶�Ӧ�������յ��������֮������ݽ��л���
%   ǰ�ĸ������ֱ��Ӧ��ʼ���x,y�������Լ��յ��x,y�����飬����Ӧ���������������ʹ��ģ�����ݣ������ʱ�����ǣ����һ����ָ����figure�ı��
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

