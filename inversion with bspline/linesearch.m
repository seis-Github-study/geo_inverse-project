function [point,sumdistance,big_num] = linesearch( row_big,row,x,y )
%����Ѱ��ĳһ���Ƿ���������㣬���û�еĻ������Ѷ�Ӧ�ķ��س�Ϊ0,ͬʱbig_num��ӳ������Ŀ��ע������
%   row_big����ĳһ�е�Ѱ�ҵ����ݣ�x��y�ֱ����������ĵ������,row������Ѱ�ҵ��е�y����
%   point������Ѱ�ҵĵ�����꣬�����������Ļ�����ֱ���x_left,y_left,x_right,y_right,���ֻ��һ���㷵�صĻ��������ǡ������ͬһ����sumdistance�������֮��
%   bignum����Ѱ�ҵ��Ĵ����Ŀ��2�����ҵ������������㣬1�����ҵ������Ϸ��������·����Ǹ��㣬0����ֻ�ҵ�����߻����ұ��Ǹ��㣬���Ǿ��������ֵ������û���ҵ���㣬��˸����޷�ʹ��
%   -1�����ҵ���һ��С����ֵ��һ���㣬ʹ�����ٽ���ֵ������Ҫע��������ַ����Ľ��ֻ����������ȫ����Ч������½���ʹ�ã�������0���ֵ����ֱ��ʹ�ñ����ٶȼ�ȥ��
%   -1�Ļ���������ڵ�һ�����λ��
%   Ϊ�˷��㣬0�ľ��뽫��Ϊ���޴�
big_num=0;
num_left=0;
num_right=0;
n=length(row_big(1,:));
point=zeros(1,4);
updist=4;%ָ���ٽ���ֵ�����ܵĸ��֮��Ĳ�࣬�������֮�����Ŀ������������֮��ľ��뼴Ϊ���룬��Ҫע�⣡����������
if abs(row-y)<10e-10 %���ж��ڲ���ͬһ��
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
else        %����ͬһ��
    if row_big(1,x)>0   %���Ƿ������ȿ��Ƿ���ͬһ����û�д��
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
     if num_left+num_right>1.5 %����������
         big_num=2;
         point(1,1:2)=big_left;
         point(1,3:4)=big_right;
     else
         if num_left>0
             if sqrt((big_left(1,1)-x)^2+(big_left(1,2)-y)^2)<=updist %��������С���ٽ���
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

