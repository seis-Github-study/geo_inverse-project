G_sensitivity=[0,2,3;0,4,4;1,2,3]
I=find(G_sensitivity<10e-15);
G_sensitivity(I)=-1;
%找到所有小点，设为0
II=find(G_sensitivity.*(G_sensitivity-2)<=10e-15);
 G_sensitivity(II)=0;