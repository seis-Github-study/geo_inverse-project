startpointx=[1,2,3,4,5,6];
startpointy=[1,2,3,4,5,6];
endpointx=[10,11,12,13];
endpointy=[10,11,12,13];
% ����ֱ�洢����ʼ���x�����y���꣬����ֻ�Ǹ����ӡ�

%�洢���ݸ���
startnum=length(startpointx);
endnum=length(endpointx);

%�����͵������ֱ����������ֵĴ�����ȷ����Χ����ȻҲ�����Լ�������
xmin=floor(min(min(startpointx),min(endpointx)));
xmax=ceil(max(max(startpointx),max(endpointx)));

ymin=floor(min(min(startpointy),min(endpointy)));
ymax=ceil(max(max(startpointy),max(endpointy)));
square=1 ;%����֮��ľ���

xnum=(xmax-xmin)/square;      %���ݷ�Χ
ynum=(ymax-ymin)/square;
%G����Ϊ��Ӧ�ľ����С��G_raypath������·���ܶȵ��Ǹ�����
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

%% ׼��������ɣ�������м���
%���ڻ��������ۻ���Ҫ������ɣ�����ⲿ��ֻ��һ��������ʾ���������������յĽ��
%�ⲿ�ֵ�Ŀ�ļ��ѽ�С�����жȾ������ֵ���ò�ֵ�������������Ҫ�ú��о�һ��˫���Բ�ֵ����㷨
 G_sensitivity(1,i)=zeros(1,xum*ynum);
for i=1:xnum*ynum
    G_sensitivity(1,i)=sum(G(:,i));
end

G_sensitivity_matrix=reshape(G_sensitivity,ynum,xnum);


%% �������G��M��������ؽ��󣬼���m��ֵ




%% ��������µ�������