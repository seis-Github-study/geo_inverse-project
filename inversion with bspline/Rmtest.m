I=find(Rm1>0.01);
nn=length(I);
for i=1:nn
Rm_new(i,1)=Rm(I(i),1);
Rm1_new(i,1)=Rm1(I(i),1);
end
figure(1)
%subplot(2,1,1)
plot(m_nonoise,'r-')
%subplot(2,1,2)
hold on
plot(m_nonoise1,'b-')
hold on
plot(m,'g--')
tt=Rm1-Rm;
figure(2)
plot(tt)