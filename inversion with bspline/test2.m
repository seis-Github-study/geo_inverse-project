G1=[0,0,0,  0,0,0,       0,0,0       0,0,0   56,0.085,12,   0,5,54,     1,5,1,     4,5,6];
G2=[0,0,0,  0,0,0,       0,0,0       0,0,0   56,0.085,12,   0,0,54,     1,5,1,     4,5,6];
G3=[0,0,0,  0,0,0,       0,0,0       0,0,0   56,0.085,12,   0,0,54,     1,5,1,     4,5,6];
G4=[0,0,0,  0,0,0,       0,0,0       0,0,0   56,0.085,12,   0,0,54,     1,5,1,     4,5,6];
G5=[0,0,0,  0,0,0,       0,0,0       0,0,0   56,0.085,12,   0,0,54,     1,5,1,     4,5,6];
G=[G1;G2;G3;G4;G5];
xnum=3;
ynum=8;
%% ׼��������ɣ�������м���
%���ڻ��������ۻ���Ҫ������ɣ�����ⲿ��ֻ��һ��������ʾ���������������յĽ��
%�ⲿ�ֵ�Ŀ�ļ��ѽ�С�����жȾ������ֵ���ò�ֵ�������������Ҫ�ú��о�һ��˫���Բ�ֵ����㷨
G_sensitivity=zeros(1,xnum*ynum);
for i=1:xnum*ynum
    G_sensitivity(1,i)=sum(G(:,i));
end


%���ж�����ֵ
G_sensitivity_low=2;
%��Ϊ0�ľ���Ԫ����Ϊ-1
I=find(G_sensitivity<10e-15);
G_sensitivity(I)=-1;
%�ҵ�����С�㣬��Ϊ0
[row,col]=find((G_sensitivity.*(G_sensitivity-G_sensitivity_low))<=10e-15);
G_sensitivity(row,col)=0;
G_sensitivity_matrix=reshape(G_sensitivity,xnum,ynum)';
clear G_sensitivity
G_sensitivity=G_sensitivity_matrix;
[row,col]=find((G_sensitivity_matrix.*(G_sensitivity_matrix-G_sensitivity_low))<=10e-15);
 %% ���м�������ֵ
 small_num=length(row);
% point,pointa,,pointb,,pointc,,pointd,wa,wb,wc,wd     (x,y)
W_small=zeros(small_num,14);
%��ֵ�ľ���֮������
updistance=10;
for i=1:1
    W_small(i,1)=col(i);
    W_small(i,2)=row(i);% ��¼��ȡ��С�������

    %a����
    %datapointline �洢x_left,y_left,x_right,y_right,distance��bignum
    datapointline=zeros(ynum,6);
    %[datapointline(row(i),1:4),distance,bignum]=linesearch(G_sensitivity(row(i),:),row(i),row(i),col(i));
    [datapointline(row(i),1:4),datapointline(row(i),5),datapointline(row(i),6)]=linesearch(G_sensitivity(row(i),:),row(i),row(i),col(i));
    %b ����
    datapointcol=zeros(xnum,6);
    [datapointcol(col(i),1:4),datapointcol(col(i),5),datapointcol(col(i),6)]=colsearch(G_sensitivity(:,col(i)),col(i),row(i),col(i));
%     %c����
%     %c1 ���õ�����������
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
%     %c2 ֻ��һ��2��������
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
%     % c3 ֻ���ٽ�������
%     %c31 ���������ٽ������
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
%     %c31 ֻ��һ���ٽ�������
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
%     % e����
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
%     %ɸѡ��Сֵ��
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
%     %ɸѡ��Сֵ��
%     if line_low>10e-15
%         [a,b]=min(datapointline(row(i)+1:ynum,5));
%         datapointlinelow=datapointline(b,:);
%     end
%     %e3
%     % f����
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
%     %ɸѡ��Сֵ��
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
%     %ѡȡ��Сֵ��
%     if col_right>10e-15
%         [a,b]=min(datapointlcol(col(i)+1:xnum,5));
%         datapointcolleft=datapointcol(b,:);
%     end
%     
%     %��Ҫע����ǣ������Ҫ���Ͼ������ֵ�жϵĻ����������line_up֮��Ĳ�����Ϊ0����
%     %��ʱ���Կ��ǲ�����Ϊ�ҹ��Ʋ���������ô����ĵ�·�߷ֲ�
%     %���������Ҫ��Ҫ���ע�⣡������������������������������������������������������������������������������������������ֵ
%     
%     % h����
%     %h1
%     if line_up*line_low*col_left*col_right>0
%         if datapointlineup(1,5)+datapointlinelow(1,5)<=datapointcolleft(1,5)+datapointright(1,5)
%             [w1,w2,w3,w4]=lineinterp4(datapointlineup(1,1:4),datapointlinelow(1,1:4),col(i),row(i));%ע���������������Ҫ�ж�һ���ĸ��ڶ���Ԫ��Ϊ0��Ȼ���������������Ȩ�ض���0
%             W_small(i,3:10)=[datapointlineup(1,1:4),datapointlinelow(1,1:4)];
%             W_small(i,11:14)=[w1,w2,w3,w4];
%             continue
%         else
%             [w1,w2,w3,w4]=colinterp4(datapointcolleft(1,1:4),datapointcolright(1,1:4),col(i),row(i));%ע���������������Ҫ�ж�һ���ĸ��ڶ���Ԫ��Ϊ0��Ȼ���������������Ȩ�ض���0
%             W_small(i,3:10)=[datapointcolleft(1,1:4),datapointcolright(1,1:4)];
%             W_small(i,11:14)=[w1,w2,w3,w4];
%             continue
%         end
%     end
%     %h2
%     if line_up*line_low>0 &&col_left*col_right<10e-16
%         [w1,w2,w3,w4]=lineinterp4(datapointlineup(1,1:4),datapointlinelow(1,1:4),col(i),row(i));%ע���������������Ҫ�ж�һ���ĸ��ڶ���Ԫ��Ϊ0��Ȼ���������������Ȩ�ض���0
%         W_small(i,3:10)=[datapointlineup(1,1:4),datapointlinelow(1,1:4)];
%         W_small(i,11:14)=[w1,w2,w3,w4];
%         continue
%     end
%     %h3
%     if line_up*line_low<10e-16 &&col_left*col_right>0
%         [w1,w2,w3,w4]=colinterp4(datapointcolleft(1,1:4),datapointcolright(1,1:4),col(i),row(i));%ע���������������Ҫ�ж�һ���ĸ��ڶ���Ԫ��Ϊ0��Ȼ���������������Ȩ�ض���0
%         W_small(i,3:10)=[datapointcolleft(1,1:4),datapointcolright(1,1:4)];
%         W_small(i,11:14)=[w1,w2,w3,w4];
%         continue
%     end
%     %h4
%     if line_up+line_low+col_left+col_right>0
%         h4data=[datapointlineup;datapointlinelow;datapointcolleft;datapointcolright];
%         [a,b]=min(hedata(:,5));
%         h4updistance=5;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!����ֵ���ǳ���Ҫ
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
%     %h5 ���������������������֮��Ĳ�࣬ѡ����С���Ǹ������������С���ҷ������ǵ�Ҫ�󣬾��������ٽ���ֵ
%     %\h6
%     %û���κξٶ������������Ȩ�ؾ����ﶼ��0��������G������Ҫ�����жϣ���������Ļ�������������������������������������������������
%     %�����ñ����ٶȽ���ȥ��
%     
 end