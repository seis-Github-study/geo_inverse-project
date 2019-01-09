function [point,sumdistance,big_num] = colsearch( col_big,col,y,x )
%����Ѱ��ĳһ���Ƿ���������㣬���û�еĻ������Ѷ�Ӧ�ķ��س�Ϊ0,ͬʱbig_num��ӳ������Ŀ��ע������
%   row_big����ĳһ�е�Ѱ�ҵ����ݣ�x��y�ֱ����������ĵ������,row������Ѱ�ҵ��е�y����
%   point������Ѱ�ҵĵ�����꣬�����������Ļ�����ֱ���x_left,y_left,x_right,y_right,���ֻ��һ���㷵�صĻ��������ǡ������ͬһ����sumdistance�������֮��
%   bignum����Ѱ�ҵ��Ĵ����Ŀ��2�����ҵ������������㣬1�����ҵ������Ϸ��������·����Ǹ��㣬0����ֻ�ҵ�����߻����ұ��Ǹ��㣬���Ǿ��������ֵ������û���ҵ���㣬��˸����޷�ʹ��
%   -1�����ҵ���һ��С����ֵ��һ���㣬ʹ�����ٽ���ֵ������Ҫע��������ַ����Ľ��ֻ����������ȫ����Ч������½���ʹ�ã�������0���ֵ����ֱ��ʹ�ñ����ٶȼ�ȥ��
%   -1�Ļ���������ڵ�һ�����λ��
%   Ϊ�˷��㣬0�ľ��뽫��Ϊ���޴�
big_num=0;
num_down=0;
num_up=0;
n=length(col_big(:,1));
point=zeros(1,4);
updist=4;%ָ���ٽ���ֵ�����ܵĸ��֮��Ĳ�࣬�������֮�����Ŀ������������֮��ľ��뼴Ϊ���룬��Ҫע�⣡����������
if abs(col-x)<10e-10 %���ж��ڲ���ͬһ��
    if y-1>0
    for j=y-1:-1:1
        if col_big(j,1)>0
            num_down=num_down+1
            big_down(1,1)=col;
            big_down(1,2)=j;
            break
        end
    end
    end
    if y+1<=n
    for j=y+1:n
        if col_big(j,1)>0
            num_up=num_up+1
            big_up(1,1)=col;
            big_up(1,2)=j;
            break
        end
    end
    end
else        %����ͬһ��
    if col_big(y,1)>0   %���Ƿ������ȿ��Ƿ���ͬһ����û�д��
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
        if y+1<=n
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
     if num_down+num_up>1.5 %����������
         big_num=2;
         point(1,1:2)=big_down;
         point(1,3:4)=big_up;
     else
         if num_down>0
             if sqrt((big_down(1,1)-x)^2+(big_down(1,2)-y)^2)<=updist %��������С���ٽ���
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
end

