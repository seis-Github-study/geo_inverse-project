function [ arg ] = R( m_data )
%计算曲率，需要注意第一列为m,第二列为GM-d
%   其实没啥影响.....
n=length(m_data(:,1));
Y=m_data(:,1);
X=m_data(:,2);
arg=zeros(n-2,2);
% arg(1,:)=nan;
% arg(2,:)=nan;
R=zeros(n,1);
for i=2:n-1
    a=sqrt((Y(i,1)-Y(i-1,1))^2+(X(i,1)-X(i-1,1))^2);
    b=sqrt((Y(i,1)-Y(i+1,1))^2+(X(i,1)-X(i+1,1))^2);
    c=sqrt((Y(i+1,1)-Y(i-1,1))^2+(X(i+1,1)-X(i-1,1))^2);
    p=(a+b+c)/2;
    R(i,1)=a*b*c/(4*sqrt(p*(p-a)*(p-b)*(p-c)));
end
[arg(:,1),arg(:,2)]=sort(R(2:n-1,1));

