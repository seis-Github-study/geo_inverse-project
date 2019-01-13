% for i=1:11
%     for j=1:10
%         m(10*(i-1)+j,1)=(-1)^(i+j);
%     end
% end
% m=random('norm',1,0.1,110,1);
for i=1:11
    for j=1:10
        m(10*(i-1)+j,1)=3+(-1)^(i+j);
    end
end
d=G*m;


%% Éú³ÉW¾ØÕó
[row_small,col_small]=find((G_sensitivity.*(G_sensitivity-G_sensitivity_low))<=10e-15);
[row_big,col_big]=find(((G_sensitivity-G_sensitivity_low))>10e-15);
nn=length(row_small);
W=eye(xnum*ynum,xnum*ynum);
for k=1:small_num
    if sum(W_small(k,11:14))<10e-10
        %d_new=d_new-ref_c*G(:,W_small(k,1)+(W_small(k,2)-1)*xnum);
        continue
    end
    for t=1:4
        if (W_small(k,2+2*t))>0
    W(W_small(k,1)+(W_small(k,2)-1)*xnum,W_small(k,1+2*t)+(W_small(k,2+2*t)-1)*xnum)=W_small(k,10+t);
    W(W_small(k,1)+(W_small(k,2)-1)*xnum,W_small(k,1)+(W_small(k,2)-1)*xnum)=0;
        end
    end
end
%%
[U,S,V]=svd(G);
Vp=V(:,1:10);
Up=U(:,1:10);
for i=1:10
    S_inverse=1/S(i,i);
end
Rm=diag(Vp*Vp');
m_nonoise=Vp*S_inverse*Up'*d;
II=find(((G_sensitivity-G_sensitivity_low))>10e-15);
diff=norm((m(II)-m_nonoise(II)));
[U1,S1,V1]=svd(G_new);
Vp1=V1(:,1:10);
Up1=U1(:,1:10);
Rm1=diag(Vp1*Vp1');
for i=1:10
    S_inverse1=1/S1(i,i);
end
m_nonoise1=Vp1*S_inverse1*Up1'*d;
diff1=sum(abs(m(II)-m_nonoise1(II)));