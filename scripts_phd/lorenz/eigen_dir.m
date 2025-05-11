sigma=10;
b=8/3;
r=24;

c=sqrt(b*(r-1));
x=c;
y=c;
z=r-1;

% x=0;
% y=0;
% z=0;

J=[-sigma, sigma,0;r-z,-1,-x;y,x,-b];

[v,u]=eig(J);
u=diag(u);

u
v
evc=real(v(:,2));
% alpha=5;
% a0=5; a1=6;
% alpha=(a0+a1)/2;
alpha=   5.906250000000000
% 
a=alpha
x=x+alpha*evc(1);
y=y+alpha*evc(2);
z=z+alpha*evc(3);
%% ti
lorenz = @(t,x) [sigma*(x(2)-x(1));r*x(1)-1*x(1)*x(3)-x(2);x(1)*x(2)-(b)*x(3)];    % Anonymous Function
% [T,X] = ode45(lorenz, [0 1.5586], [-6.2262, -11.0027,13.0515]);
[T,X] = ode45(lorenz, [0 100], [x,y,z]);
close all; plot(T,X); grid on; grid minor;

%% analysis
z=fft(X(:,3));
f=0:1/T(end):1/mean(diff(T));
% close all; plot(f(2:end),abs(z(2:end))); grid on; grid minor; 
close all; semilogy(f(2:end),abs(z(2:end))); grid on; grid minor; 
[a,d]=max(abs(z(2:end)));
per=1/f(d+1)
%%
close all; plot3(X(end-round(per/mean(diff(T))*2):end,1),X(end-round(per/mean(diff(T))*2):end,2),X(end-round(per/mean(diff(T))*2):end,3))
np=per/mean(diff(T))
%% gen of upo data
% np=20;
np=np+1
lorenz = @(t,x) [sigma*(x(2)-x(1));r*x(1)-1*x(1)*x(3)-x(2);x(1)*x(2)-(b)*x(3)];    % Anonymous Function
% [T,X] = ode45(lorenz, [0 1.5586], [-6.2262, -11.0027,13.0515]);
[T,X] = ode45(lorenz, [0:per/np:per], [X(end,1),X(end,2),X(end,3)]);
close all; plot(mod(T,per),X(:,1:2)); 
%%
% U=[];
% R=[];
% 
% for r=0.1:0.1:4
% J=[-sigma, sigma,0;r-z,-1,-x;y,x,-b];
% 
% [v,u]=eig(J);
% u=diag(u);
% 
% U=[U,u];
% R=[R,r];
% end
% 
% plot(R,U')