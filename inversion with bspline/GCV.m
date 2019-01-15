function [ g ] =GCV( G,m,d,L,a)
%使用GCV来计算对应的值
%   G代表正演矩阵，data代表数据，L代表对应的模型正则化矩阵，a为阿尔法的数值
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

