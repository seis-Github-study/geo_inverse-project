function [ w1,w2 ] = colinterp2( point,x,y )
%lineinterp2 一维行行方向的线性插值
%   point x_left,y_left,x_right,y_right
x1=point(1,1);
y1=point(1,2);
x2=point(1,3);
y2=point(1,4);
w1=(y2-y)/(y2-y1);
w2=(y-y1)/(y2-y1);

