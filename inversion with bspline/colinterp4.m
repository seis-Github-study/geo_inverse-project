function [ w1,w2,w3,w4 ] = colinterp4( pointa,pointb,x,y )
%由行搜索得到的四组点插值得到的点
%   此处显示详细说明
x1=pointa(1,1);
y1=pointa(1,2);
x2=pointa(1,3);
y2=pointa(1,4);

x3=pointb(1,1);
y3=pointb(1,2);
x4=pointb(1,3);
y4=pointb(1,4);

if x1*x2*x3*x4>0
    w1=((x3-x)/(x3-x1))*(y2-y)/(y2-y1);
    w2=((x3-x)/(x3-x1))*(y-y1)/(y2-y1);
    w3=((x-x1)/(x3-x1))*(y4-y)/(y4-y3);
    w4=((x-x1)/(x3-x1))*(y-y3)/(y4-y3);
else
    if x1*y1==0
        error('the function of point is error')
    end
    if x2*y2==0
        w1=((y3-y)/(y3-y1));
        w2=0;
        w3=((y-y1)/(y3-y1))*(x4-x)/(x4-x3);
        w4=((y-y1)/(y3-y1))*(x-x3)/(x4-x3);
    end
    if x3*y3==0
        error('the function of point is error')
    end    
    if x4*y4==0
    w1=((y3-y)/(y3-y1))*(x2-x)/(x2-x1);
    w2=((y3-y)/(y3-y1))*(x-x1)/(x2-x1);
    w3=((y-y1)/(y3-y1));
    w4=0;
    end
end
