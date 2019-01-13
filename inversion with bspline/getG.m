function [ G,G_raypath ] = getG( xpoint,ypoint,square, xmin,xmax,ymin,ymax)
%����G����
%   ��������û�����⣬����Ҫ��һ���Ĳ���
% clear
% xpoint=[1,2];%start point and end point
% ypoint=[3,1];

% square=0.5; %���ֵ������С

% xmin=0;        %��γ�ȵķ�Χ
% xmax=3;
% ymin=0;
% ymax=3;

xnum=(xmax-xmin)/square;      %���ݷ�Χ
ynum=(ymax-ymin)/square;
xline=xmin:square:xmax;
yline=ymin:square:ymax; %�ߵķ�Χ

G_raypath=zeros(ynum,xnum);
G_raylength=zeros(ynum,xnum);
% G���������Ϊx��������Ϊy����x=j,y=i;
if (xpoint(2)-xpoint(1))~=0
k=(ypoint(2)-ypoint(1))/(xpoint(2)-xpoint(1)); %����б��
else
    %k is not exsit
end

% if k==0
%     %k is zero
% end
if k>=0
    if ypoint(2)<ypoint(1)          %if the start point is in the higher value,we change the start and end point which has different effect on the result
        x_tmp=xpoint(1);
        y_tmp=ypoint(1);
        xpoint(1)=xpoint(2);
        ypoint(1)=ypoint(2);
        xpoint(2)=x_tmp;
        ypoint(2)=y_tmp;
    end
    x1=floor(xpoint(1)/square)*square;%���㿪ʼ��������һ������
    y1=floor(ypoint(1)/square)*square;
    x_start=find(xline==x1);
    y_start=find(yline==y1);
    
    x2=floor(xpoint(2)/square)*square;%���������������һ������
    y2=floor(ypoint(2)/square)*square;
    if (xpoint(2)/square)-floor(xpoint(2)/square)==0
        x_end=find(xline==x2)-1;
    else
        x_end=find(xline==x2);
    end
    if (ypoint(2)/square)-floor(ypoint(2)/square)==0
        y_end=find(yline==y2)-1;
    else
        y_end=find(yline==y2);
    end
    % y_end=find(yline==y2)-1;
    
    x_now=x_start;
    y_now=y_start;
    
    x1=xpoint(1); %start location of the point(x,y)
    y1=ypoint(1);
    
    while(~((x_now==x_end) && (y_now==y_end)))
        x_next_grid=xline(x_now)+square;
        y_next_grid=yline(y_now)+square;
        k_tmp=(y_next_grid-y1)/(x_next_grid-x1);
        if k_tmp>k
            x_next=x_now+1;       % location of the next grid
            y_next=y_now;
            x2=xline(x_next);
            y2=y1+k*(x2-x1);
        else if k_tmp==k
                x_next=x_now+1;       % location of the next grid
                y_next=y_now+1;
                x2=xline(x_next);
                y2=y1+k*(x2-x1);
            else
                y_next=y_now+1;
                x_next=x_now;
                y2=yline(y_next);
                x2=(y2-y1)/k+x1;
            end
        end
        G_raypath(y_now,x_now)=1 ;%whether go through the grid
        G_raylength(y_now,x_now)=sqrt((y2-y1)^2+(x2-x1)^2);%G sensivity matrix
        
        x1=x2;% update information
        y1=y2;
        x_now=x_next;
        y_now=y_next;
    end
    G_raypath(y_next,x_next)=1;%whether go through the grid
    G_raylength(y_next,x_next)=sqrt((y2-ypoint(2))^2+(x2-xpoint(2))^2);%G sensivity matrix
    
    x1=x2;% update information
    y1=y2;
    x_now=x_next;
    y_now=y_next;
else                                % k<0
    if ypoint(2)>ypoint(1)          %if the start point is in the higher value,we change the start and end point which has different effect on the result
        x_tmp=xpoint(1);
        y_tmp=ypoint(1);
        xpoint(1)=xpoint(2);
        ypoint(1)=ypoint(2);
        xpoint(2)=x_tmp;
        ypoint(2)=y_tmp;
    end
    x1=floor(xpoint(1)/square)*square;%���㿪ʼ��������һ������
    y1=floor(ypoint(1)/square)*square;
    
    
    x2=floor(xpoint(2)/square)*square;%���������������һ������
    y2=floor(ypoint(2)/square)*square;
    x_start=find(xline==x1);
    y_end=find(yline==y2);
    if (xpoint(2)/square)-floor(xpoint(2)/square)==0
        x_end=find(xline==x2)-1;
    else
        x_end=find(xline==x2);
    end
    if (ypoint(1)/square)-floor(ypoint(1)/square)==0
        y_start=find(yline==y1)-1;
    else
        y_start=find(yline==y1);
    end
    
    x_now=x_start;
    y_now=y_start;
    
    x1=xpoint(1); %start location of the point(x,y)
    y1=ypoint(1);
    
    while(~(x_now==x_end &&y_now==y_end))
        x_next_grid=xline(x_now)+square;
        y_next_grid=yline(y_now);
        k_tmp=(y_next_grid-y1)/(x_next_grid-x1);
        if k_tmp<k
            x_next=x_now+1;       % location of the next grid
            y_next=y_now;
            x2=xline(x_next);
            y2=y1+k*(x2-x1);
        else if  k_tmp==k
                x_next=x_now+1;       % location of the next grid
                y_next=y_now-1;
                x2=xline(x_next);
                y2=y1+k*(x2-x1);
            else
                y_next=y_now-1;
                x_next=x_now;
                y2=yline(y_next)+square;
                x2=(y2-y1)/k+x1;
            end
        end
        G_raypath(y_now,x_now)=1; %whether go through the grid
        G_raylength(y_now,x_now)=sqrt((y2-y1)^2+(x2-x1)^2);%G sensivity matrix
        
        x1=x2;% update information
        y1=y2;
        x_now=x_next;
        y_now=y_next;
    end
    G_raypath(y_next,x_next)=1; %whether go through the grid
    G_raylength(y_next,x_next)=sqrt((y2-ypoint(2))^2+(x2-xpoint(2))^2);%G sensivity matrix
    
    x1=x2;% update information
    y1=y2;
    x_now=x_next;
    y_now=y_next;
end
  G(1,:)=reshape(G_raylength',1,xnum*ynum);
  
   

end

