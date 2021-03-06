clear
col_big=[0,0,0,  4,4,0,  1,2,3,  0,0,0,  0,0,0,  0,0]';
col=3;
x=5;
y=6;

big_num=0;
num_down=0;
num_up=0;
n=length(col_big(:,1));
point=zeros(1,4);
updist=4;%指最临近插值所接受的格点之间的差距，即网格点之间的数目，乘以网格点的之间的距离即为距离，需要注意！！！！！！
if abs(col-x)<10e-10 %先判断在不在同一列
    if y-1>0
    for j=y-1:-1:1
        if col_big(j,1)>0
            num_down=num_down+1;
            big_down(1,1)=col;
            big_down(1,2)=j;
            break
        end
    end
    end
    if y+1<n
    for j=y+1:n
        if col_big(j,1)>0
            num_up=num_up+1;
            big_up(1,1)=col;
            big_up(1,2)=j;
            break
        end
    end
    end
else        %不在同一列
    if col_big(y,1)>0   %看是否在首先看是否在同一列有没有大点
        big_num=1;
        point(1,1)=col;
        point(1,2)=y;
        point(1,3:4)=0;
    else
        if y-1>0
        for j=y-1:-1:1
            if col_big(j,1)>0
                num_down=num_down+1;
                big_down(1,2)=j;
                big_down(1,1)=col;
                break
            end
        end
        end
        if y+1<n
        for j=y+1:n
            if col_big(j,1)>0
                num_up=num_up+1;
                big_up(1,2)=j;
                big_up(1,1)=col;
                break
            end
        end
        end
    end   
end
     if num_down+num_up>1.5 %即有两个点
         big_num=2;
         point(1,1:2)=big_down;
         point(1,3:4)=big_up;
     else
         if num_down>0
             if sqrt((big_down(1,1)-x)^2+(big_down(1,2)-y)^2)<=updist %即距离最小的临近点
                 point(1,1:2)=big_down;
                 point(1,3:4)=0;
                 big_num=-1;
             end
         end
         if num_up>0
             if sqrt((big_up(1,1)-x)^2+(big_up(1,2)-y)^2)<=updist
                 point(1,1:2)=big_up;
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

