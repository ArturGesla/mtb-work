clc; clear;
x=-2:1:2; y=[-3:1:3]'; 


z=x.^2-(y-3).^2;  zh=@(x,y) x^2-(y-3)^2;

x=-2:1/10:2; y=[-3:1/10:3]'; z=1/2*x.^2+1/4*y.^4-1/2*y.^2; zh=@(x,y) 1/2*x.^2+1/4*y.^4-1/2*y.^2;
contour(x,y,z);
% surf(x,y,z);

%%

clc;
x=[5;5];
% x=[0.2;0.2];
xa=[x];
eps=1e-6;


%%

d2fdx2=(zh(x(1)+eps,x(2))-2*zh(x(1),x(2))+zh(x(1)-eps,x(2)))/eps^2;
d2fdxdy=(zh(x(1)+eps,x(2)+eps)+zh(x(1)-eps,x(2)-eps)-zh(x(1)+eps,x(2)-eps)-zh(x(1)-eps,x(2)+eps))/eps^2/4;
d2fdy2=(zh(x(1),x(2)+eps)-2*zh(x(1),x(2))+zh(x(1),x(2)-eps))/eps^2;

jac=[d2fdx2,d2fdxdy; d2fdxdy, d2fdy2];
g=[(zh(x(1)+eps,x(2))-zh(x(1)-eps,x(2)))/2/eps;(zh(x(1),x(2)+eps)-zh(x(1),x(2)-eps))/2/eps];
x=x-jac\g;
xa=[xa,x];
det(jac)
norm(g)
%%
hold on;
plot(xa(1,:),xa(2,:),'-x')