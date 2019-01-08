function [point,sumdistance,big_num] = linesearch( row_big,row,x,y )
%用来寻找某一行是否有两个大点，如果没有的话，将把对应的返回成为0,同时big_num反映大点的数目，注意区分
%   row_big代表某一行的寻找的数据，x，y分别代表所计算的点的坐标,row代表所寻找的行的y坐标
%   point代表所寻找的点的坐标，如果是两个点的话，则分别是x_left,y_left,x_right,y_right,如果只有一个点返回的话，则代表恰好是在同一列上sumdistance代表距离之和
%   bignum代表寻找到的大点数目，2代表找到了左右两个点，1代表找到了正上方或者正下方的那个点，0代表只找到了左边或者右边那个点，但是距离大于阈值，或者没有找到大点，因此根本无法使用
%   -1代表找到了一个小于阈值的一个点，使用最临近插值，但需要注意的是这种方法的结果只在其他方法全部无效的情况下进行使用，而对于0这个值，则直接使用背景速度减去。
%   -1的话这个点会放在第一个点的位置
%   为了方便，0的距离将设为无限大
big_num=0;
num_left=0;
num_right=0;
n=length(row_big(1,:));
point=zeros(1,4);
updist=4;%指最临近插值所接受的格点之间的差距，即网格点之间的数目，乘以网格点的之间的距离即为距离，需要注意！！！！！！
if abs(row-y)<10e-10 %先判断在不在同一行
    if x-1>0
    for j=x-1:-1:1
        if row_big(1,j)>0
            num_left=num_left+1;
            big_left(1,1)=j;
            big_left(1,2)=row;
            break
        end
    end
    end
    if x+1<n
    for j=x+1:n
        if row_big(1,j)>0
            num_right=num_right+1;
            big_right(1,1)=j;
            big_right(1,2)=row;
            break
        end
    end
    end
else        %不在同一行
    if row_big(1,x)>0   %看是否在首先看是否在同一列有没有大点
        big_num=1;
        point(1,1)=x;
        point(1,2)=row;
        point(1,3:4)=0;
    else
        if x-1>0
        for j=x-1:-1:1
            if row_big(1,j)>0
                num_left=num_left+1;
                big_left(1,1)=j;
                big_left(1,2)=row;
                break
            end
        end
        end
        if x+1<n
        for j=x+1:n
            if row_big(1,j)>0
                num_right=num_right+1;
                big_right(1,1)=j;
                big_right(1,2)=row;
                break
            end
        end
        end
    end   
end
     if num_left+num_right>1.5 %即有两个点
         big_num=2;
         point(1,1:2)=big_left;
         point(1,3:4)=big_right;
     else
         if num_left>0
             if sqrt((big_left(1,1)-x)^2+(big_left(1,2)-y)^2)<=updist %即距离最小的临近点
                 point(1,1:2)=big_left;
                 point(1,3:4)=0;
                 big_num=-1;
             end
         end
         if num_right>0
             if sqrt((big_right(1,1)-x)^2+(big_right(1,2)-y)^2)<=updist
                 point(1,1:2)=big_right;
                 point(1,3:4)=0;                 
                 big_num=-1;
             end
         end
     end
     if big_num==2
         sumdistance=sqrt((point(1,1)-x)^2+(point(1,2)-y)^2)+sqrt((point(1,3)-x)^2+(point(1,4)-y)^2);
     end
     if abs(big_num)==1
         sumdistance=sqrt((point(1,1)-x)^2+(point(1,2)-y)^2);
     end
     if big_num==0
         sumdistance=inf;
     end
end

