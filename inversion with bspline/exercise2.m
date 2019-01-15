%% problem a 
clear
load('ifk.mat')
x=0.025:0.05:0.975;
x=x';
y=x; %生成数据
G=zeros(20,20);
for i=1:20
    for j=1:20
        G(i,j)=x(j,1)*exp(-x(j,1)*y(i,1))*0.05;
    end
end

%% problem b
misfit=sqrt(20)*10^-4;

%% preparation
[U,S,V]=svd(G);
%rank(G)=8 从之前的题目结论可以看到
% S_inverse=zeros(20,20);
% for i=1:8
%     S_inverse(i,i)=1/S(i,i);
% end

%% problem c
L0=zeros(20,20);
for i=1:20
    L0(i,i)=1;
end
% 先得到一个初始α，然后乘以一定的量级再算
for i=1:11
    alpha_try(i,1)=10^(i-6); %-5到5次方
    m=getm(G,d,L0,alpha_try(i,1));
    alpha_try(i,2)=norm(L0*m);
    alpha_try(i,3)=norm(G*m-d);
end
%结果应该在10-5到100之间
% L-CURVE
number=5000;
low=-6;
up=2;
interval=(up-low)/number;
t=low:interval:up;
%a_c_L=zeros(1000,1);
m_c_L0=zeros(number,3);
m_L0_value=zeros(20,number);
for i=1:number
    m_c_L0(i,1)=10^t(i);%α的数值
    m_L0_value(:,i)=getm(G,d,L0,m_c_L0(i,1));
    m_c_L0(i,2)=norm(L0*m_L0_value(:,i));%||m||
    m_c_L0(i,3)=norm(G* m_L0_value(:,i)-d); %||Gm-d||
end
figure(1)
plot(log10(m_c_L0(:,3)),(m_c_L0(:,2)))
xlabel('||Gm-d||')
ylabel('||m||')
figure(2)
plot(m_c_L0(:,3),m_c_L0(:,2))
xlabel('||Gm-d||')
ylabel('||m||')
xlim([0,0.005])
%m=getm(G,d,L0,0.0000457);
%a=6.792036326171850e-04
m=getm(G,d,L0,6.792036326171850e-04);
norm(L0*m)
norm(G*m-d)
figure(3)
plot(x,m)
xlabel('x')
ylabel('m')
title('a=6.792*10^-4')
%%
% GCV
g0=zeros(number,2);
for i=1:number
g0(i,1)=m_c_L0(i,1);
g0(i,2)=GCV( G,m_L0_value(:,i),d,L0, m_c_L0(i,1));
end

figure(4)
plot(log10(g0(:,1)),log10(g0(:,2)))
[a,b]=min(g0(:,2));%得到最小值
m=getm(G,d,L0,g0(b,1));
figure(5)
plot(x,m)
xlabel('x')
ylabel('m')
title('a=5.6494*10^-5')

%%
%disperancy principle
%m_c_L0第1784行 α=7.1252e-04 misfit=4.7158e-4
m=getm(G,d,L0,7.125248247522675e-04);
figure(6)
plot(x,m)
xlabel('x')
ylabel('m')
title('a=4.7158e-4')
%% 问题d
L1=zeros(19,20);
for i=1:19
    L1(i,i)=-1;
    L1(i,i+1)=1;
end
% L-CURVE
number=5000;
low=-6;
up=2;
interval=(up-low)/number;
t=low:interval:up;
%a_c_L=zeros(1000,1);
m_c_L1=zeros(number,3);
m_L1_value=zeros(20,number);
for i=1:number
    m_c_L1(i,1)=10^t(i);%α的数值
    m_L1_value(:,i)=getm(G,d,L1,m_c_L1(i,1));
    m_c_L1(i,2)=norm(L1*m_L1_value(:,i));%||m||
    m_c_L1(i,3)=norm(G* m_L1_value(:,i)-d); %||Gm-d||
end
figure(7)
plot(log10(m_c_L1(:,3)),(m_c_L1(:,2)))
xlabel('||Gm-d||')
ylabel('||m||')
figure(8)
plot(m_c_L1(:,3),m_c_L1(:,2))
xlabel('||Gm-d||')
ylabel('||m||')
%xlim([0,0.005])
%m=getm(G,d,L1,9.752800697426161e-05);
%a=9.752800697426161e-05
m=getm(G,d,L1,9.752800697426161e-05);
% norm(L0*m)
% norm(G*m-d)
figure(9)
plot(x,m)
xlabel('x')
ylabel('m')
title('a=9.7528*10^-5')
%% GCV
g1=zeros(number,2);
for i=1:number
g1(i,1)=m_c_L1(i,1);
g1(i,2)=GCV( G,m_L1_value(:,i),d,L1, m_c_L1(i,1));
end
figure(10)
plot(log10(g1(:,1)),log10(g1(:,2)))

[a,b]=min(g1(:,2));%得到最小值
%%因为这个是单调增加，然后猛地增加，根据选取的最小值10^-6解出来非常振荡，因此不可能往左边走，同时也不可能往右边走，因此只能是选取猛然增加的那个点，需要在新版的上面进行改进。
m=getm(G,d,L1,g1(b,1));
figure(11)
plot(x,m)
xlabel('x')
ylabel('m')
title('a=2.3249e-04')

%% DISPERANCY PRINCIPLE
%disperancy principle
%m_c_L0第2462行 α=0.0168 misfit=4.4663e-4
m=getm(G,d,L0,0.016811249744770);
figure(12)
plot(x,m)
xlabel('x')
ylabel('m')
title('a=0.0168')
%%
L2=zeros(18,20);
for i=1:18
    L2(i,i)=1;
    L2(i,i+1)=-2;
    L2(i,i+2)=1;
end
% L-CURVE
number=5000;
low=-6;
up=2;
interval=(up-low)/number;
t=low:interval:up;
%a_c_L=zeros(1000,1);
m_c_L2=zeros(number,3);
m_L2_value=zeros(20,number);
for i=1:number
    m_c_L2(i,1)=10^t(i);%α的数值
    m_L2_value(:,i)=getm(G,d,L2,m_c_L2(i,1));
    m_c_L2(i,2)=norm(L2*m_L2_value(:,i));%||m||
    m_c_L2(i,3)=norm(G* m_L2_value(:,i)-d); %||Gm-d||
end
figure(13)
plot(log10(m_c_L2(:,3)),(m_c_L2(:,2)))
xlabel('||Gm-d||')
ylabel('||m||')
figure(14)
plot(m_c_L2(:,3),m_c_L2(:,2))
xlabel('||Gm-d||')
ylabel('||m||')
%xlim([0,0.005])
%a=6.852306199651597e-04
 m=getm(G,d,L2,6.852306199651597e-04);
 figure(15)
 plot(x,m)
 xlabel('x')
 ylabel('m')
 title('a=6.8523*10^-4')
 
 %%
 %% GCV
g2=zeros(number,2);
for i=1:number
g2(i,1)=m_c_L2(i,1);
g2(i,2)=GCV( G,m_L2_value(:,i),d,L2, m_c_L2(i,1));
end
figure(16)
plot(log10(g2(:,1)),log10(g2(:,2)))

[a,b]=min(g2(:,2));%得到最小值
m=getm(G,d,L2,0.002744100418257);
figure(17)
plot(x,m)
xlabel('x')
ylabel('m')
title('a=0.002744')

%% 
%disperancy principle
%m_c_L0第2462行 α=0.0168 misfit=4.4663e-4
% number=5000;
% low=5;
% up=20;
% interval=(up-low)/number;
% t=low:interval:up;
% %a_c_L=zeros(1000,1);
% m_c_L2_disp=zeros(number,3);
% m_L2_value_disp=zeros(20,number);
% for i=1:number
%     m_c_L2_disp(i,1)=10^t(i);%α的数值
%     m_L2_value_disp(:,i)=getm(G,d,L2,m_c_L2_disp(i,1));
%     m_c_L2_disp(i,2)=norm(L0*m_L2_value_disp(:,i));%||m||
%     m_c_L2_disp(i,3)=norm(G* m_L2_value_disp(:,i)-d); %||Gm-d||
% end
%%
m=getm(G,d,L2,1);%无法达到所谓的误差值，后面m的值就不变化了
figure(18)
plot(x,m)
xlabel('x')
ylabel('m')
title('a=1')


%% 分辨率矩阵

%% L0
