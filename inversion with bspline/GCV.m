function [ g ] =GCV( G,m,d,L,a)
%ʹ��GCV�������Ӧ��ֵ
%   G�������ݾ���data�������ݣ�L�����Ӧ��ģ�����򻯾���aΪ����������ֵ
%G_proinverse=(G'*G+a*a.*(L'*L))\(G');
[U,V,X,C,S]=gsvd(G,L);
%Y=inv(X');
T=C'*C+a^2.*S'*S;
m1=min(length(T(:,1)),length(T(1,:)));
for i=1:m1
    if abs(T(i,i))<10^-15
        T(i,i)=0;
    else
        T(i,i)=1/T(i,i);
    end
end
G_proinverse=X'\(T*C'*U');
GG_proinv=G*G_proinverse;
n=length(GG_proinv(:,1));
I=eye(n,n);
tr=(trace(I-GG_proinv)).^2;
g=n*(norm(G*m-d)).^2/tr;
warning('off')


end

