%%
clc; clear;
%
amp=0.3;
A=[1,0,0,0,0;
    1,-2,1,0,0;
    0,1,-2,1,0;
    0,0,1,-2,1;
    0,0,0,0,1   ];
c=[0,0,0,0,0;
    -1,0,1,0,0;
    0,-1,0,1,0;
    0,0,-1,0,1;
    0,0,0,0,0];A=A+c*amp;
a2=A(2:end-1,2:end-1);

a3=[1,1,0,0,0;
    1,-2,1,0,0;
    0,1,-2,1,0;
    0,0,1,-2,1;
    0,0,0,1,1];
b3=[0;1;1;1;0];
u3=a3\b3;

a4=[  -3,1,0;
    1,-2,1;
    0,1,-3];
b4=[1;1;1];
u4=a4\b4;
%
% eig(A)
eig(a3)
eig(a4)
a1=A;
% a11=a1-2*c*0.1;
a11=A-2*c*amp;
%%
format shortE
format shortG
[evc3,evs3]=eig(a1); 
% [evc3t,evs3t]=eig(a1');
[evc3t,evs3t]=eig(a11);

[a,b1]=sort(diag(evs3));
[a,b2]=sort(diag(evs3t));
diag(evs3(b1,b1))'
diag(evs3t(b2,b2))'
%%
% evc3*evc3t' % rothogonal
evc3(:,b1)'*evc3t(:,b2) % rothogonal
% evc3(:,b1)'*evc3(:,b1) % rothogonal