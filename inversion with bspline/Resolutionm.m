function [ Rm] = Resolutionm( G,L,a )
%�����Ӧ�ķֱ��ʾ���
%   ûɶ�ý��ܵģ��ͽ�������a�ǰ�����
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
Rm=G_proinverse*G;

end

