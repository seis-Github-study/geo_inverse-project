function [ m ] = getm( G,d,L,a )
%���������Ӧ��TR��
% G�������ݾ���data�������ݣ�L�����Ӧ��ģ�����򻯾���aΪ����������ֵ
rank=8;
[U,V,X,C,S]=gsvd(G,L);
%Y=inv(X');
n=length(C(:,1));
% for i=n-rank:1
%     S(i,i)=0;
% end
T=C'*C+a^2.*S'*S;
m=min(length(T(:,1)),length(T(1,:)));
for i=1:m
    if abs(T(i,i))<10^-15
        T(i,i)=0;
    else
        T(i,i)=1/T(i,i);
    end
end
x=T*C'*U'*d;
%m=(G'*G+a*a.*(L'*L))\(G'*d);
%m=Y*x;
m=X'\x;
%m=m1-m2;
warning('off')
end

