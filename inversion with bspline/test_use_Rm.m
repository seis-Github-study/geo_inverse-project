%������÷ֱ��ʾ�����в�ֵ
clear
x=10:10:100;
x=x-5;
G=zeros(900,110);
k=1;
for m=1:3
    for n=1:3
        for i=1:10
            for j=1:10
                ypoint=[500*(m-1)+x(i),500*(n-1)+x(j)];
                xpoint=[0,1000];
                [G1,G2]=getG(xpoint,ypoint,100,0,1000,0,1100);
                G(k,:)=G1;
                k=k+1;
            end
        end
    end
end
[U,S,V]=svd(G);
Vp=V(:,1:39);
Rm=diag(Vp*Vp');
xnum=10;
ynum=11;
%% ׼��������ɣ�������м���
%���ڻ��������ۻ���Ҫ������ɣ�����ⲿ��ֻ��һ��������ʾ���������������յĽ��
%�ⲿ�ֵ�Ŀ�ļ��ѽ�С�����жȾ������ֵ���ò�ֵ�������������Ҫ�ú��о�һ��˫���Բ�ֵ����㷨
G_sensitivity=zeros(1,xnum*ynum);
%% �ش�ı�
for i=1:xnum*ynum
    G_sensitivity(1,i)=Rm(i,1);
end
%%

%���ж�����ֵ
G_sensitivity_low=0.2;
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
nearupdistance=5; %�ٽ���ֵ������
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
    %c����
    %c1 ���õ�����������
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
    %c2 ֻ��һ��2��������
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
    % c3 ֻ���ٽ�������
    %c31 ���������ٽ������
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
    %c31 ֻ��һ���ٽ�������
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
    % e����
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
                datapointline(j,5)=datapointline(j,5)*4;  %!!!!����һ�ֿ��Ƶķ����������ٽ����ֵ������ȫ�ܵ���ӭ�����ǲ��̫Զ
            end
            if datapointline(j,6)>10e-15
                line_up=line_up+1;
            end
        end
        %ɸѡ��Сֵ��
        if line_up>10e-15
            [a,b]=min(datapointline(1:row(i)-1,5));
            datapointlineup=datapointline(b,:);
        end
    end
    %e2
    if row(i)<ynum %��ֹ�����һ��
        for j=row(i)+1:ynum
            [datapointline(j,1:4),datapointline(j,5),datapointline(j,6)]=linesearch(G_sensitivity(j,:),j,row(i),col(i));
            if abs(datapointline(j,6)-2)<10e-15
                datapointline(j,5)=datapointline(j,5)./2;
            end
            if abs(datapointline(j,6)+1)<10e-15
                datapointline(j,5)=datapointline(j,5)*4;  %!!!!����һ�ֿ��Ƶķ����������ٽ����ֵ������ȫ�ܵ���ӭ�����ǲ��̫Զ
            end
            if (datapointline(j,6))>10e-15
                line_low=line_low+1;
            end
        end
        %ɸѡ��Сֵ��
        if line_low>10e-15
            [a,b]=min(datapointline(row(i)+1:ynum,5));
            datapointlinelow=datapointline(b+row(i),:);
        end
    end
    
    %e3
    % f����
    %f1
    if col(i)>0
        for j=1:col(i)-1
            [datapointcol(j,1:4),datapointcol(j,5),datapointcol(j,6)]=colsearch(G_sensitivity(:,j),j,row(i),col(i));
            if abs(datapointcol(j,6)-2)<10e-15
                datapointcol(j,5)=datapointcol(j,5)./2;
            end
            if abs(datapointcol(j,6)+1)<10e-15
                datapointcol(j,5)=datapointcol(j,5)*4;  %!!!!����һ�ֿ��Ƶķ����������ٽ����ֵ������ȫ�ܵ���ӭ�����ǲ��̫Զ
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
    %ɸѡ��Сֵ��
    
    %f2  
    if col(i)<xnum
        for j=col(i)+1:xnum
            [datapointcol(j,1:4),datapointcol(j,5),datapointcol(j,6)]=colsearch(G_sensitivity(:,j),j,row(i),col(i));
            if abs(datapointcol(j,6)-2)<10e-15
                datapointcol(j,5)=datapointcol(j,5)./2;
            end
            if abs(datapointcol(j,6)+1)<10e-10
                datapointcol(j,5)
                datapointcol(j,5)=datapointcol(j,5)*4;  %!!!!����һ�ֿ��Ƶķ����������ٽ����ֵ������ȫ�ܵ���ӭ�����ǲ��̫Զ
            end
            if (datapointcol(j,6))>10e-15
                
                col_right=col_right+1;
            end
            %ѡȡ��Сֵ��
            if col_right>10e-15
                [a,b]=min(datapointcol(col(i)+1:xnum,5));
                datapointcolright=datapointcol(b+col(i),:);
            end
        end
    end
    %     %��Ҫע����ǣ������Ҫ���Ͼ������ֵ�жϵĻ����������line_up֮��Ĳ�����Ϊ0����
    %     %��ʱ���Կ��ǲ�����Ϊ�ҹ��Ʋ���������ô����ĵ�·�߷ֲ�
    %     %���������Ҫ��Ҫ���ע�⣡������������������������������������������������������������������������������������������ֵ
    %
    % h����
    %h1
    if line_up*line_low*col_left*col_right>0
        if datapointlineup(1,5)+datapointlinelow(1,5)<=datapointcolleft(1,5)+datapointcolright(1,5)
            [w1,w2,w3,w4]=lineinterp4(datapointlineup(1,1:4),datapointlinelow(1,1:4),col(i),row(i));%ע���������������Ҫ�ж�һ���ĸ��ڶ���Ԫ��Ϊ0��Ȼ���������������Ȩ�ض���0
            W_small(i,3:10)=[datapointlineup(1,1:4),datapointlinelow(1,1:4)];
            W_small(i,11:14)=[w1,w2,w3,w4];
            continue
        else
            [w1,w2,w3,w4]=colinterp4(datapointcolleft(1,1:4),datapointcolright(1,1:4),col(i),row(i));%ע���������������Ҫ�ж�һ���ĸ��ڶ���Ԫ��Ϊ0��Ȼ���������������Ȩ�ض���0
            W_small(i,3:10)=[datapointcolleft(1,1:4),datapointcolright(1,1:4)];
            W_small(i,11:14)=[w1,w2,w3,w4];
            continue
        end
    end
    %h2
    if line_up*line_low>0 &&col_left*col_right<10e-16
        [w1,w2,w3,w4]=lineinterp4(datapointlineup(1,1:4),datapointlinelow(1,1:4),col(i),row(i));%ע���������������Ҫ�ж�һ���ĸ��ڶ���Ԫ��Ϊ0��Ȼ���������������Ȩ�ض���0
        W_small(i,3:10)=[datapointlineup(1,1:4),datapointlinelow(1,1:4)];
        W_small(i,11:14)=[w1,w2,w3,w4];
        continue
    end
    %     %h3
    if line_up*line_low<10e-16 &&col_left*col_right>0
        [w1,w2,w3,w4]=colinterp4(datapointcolleft(1,1:4),datapointcolright(1,1:4),col(i),row(i));%ע���������������Ҫ�ж�һ���ĸ��ڶ���Ԫ��Ϊ0��Ȼ���������������Ȩ�ض���0
        W_small(i,3:10)=[datapointcolleft(1,1:4),datapointcolright(1,1:4)];
        W_small(i,11:14)=[w1,w2,w3,w4];
        continue
    end
    %h4 ������ֵ�õ��������(�л��в�ֵ�õ�)
    if line_up+line_low+col_left+col_right>0
        h4data=[datapointlineup;datapointlinelow;datapointcolleft;datapointcolright];
        [a,b]=min(hedata(:,5));
        h4updistance=5;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!����ֵ���ǳ���Ҫ
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
    %     ���������������������֮��Ĳ�࣬ѡ����С���Ǹ������������С���ҷ������ǵ�Ҫ�󣬾��������ٽ���ֵ,��Ϊ��linesearch��colsearch�����Ѿ��ӹ������ˣ����������Ͳ��������������Ѿ��ӹ�
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
    h2updistance=5;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!����ֵ���ǳ���Ҫ
    if a>updistance
        break
    end
    W_small(i,3:4)=h2data(b,1:2);
    W_small(i,11)=1;
    
    
end

%     %\h6
%     %û���κξٶ������������Ȩ�ؾ����ﶼ��0��������G������Ҫ�����жϣ���������Ļ�������������������������������������������������
%     %�����ñ����ٶȽ���ȥ��
%% ����Ȩ�ؾ���W_small�������飨ע����ֻ�Ǽ�¼����Ϣ���������¾���w��
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
[U1,S1,V1]=svd(G_new);
Vp1=V1(:,1:39);
Rm1=diag(Vp1*Vp1');
