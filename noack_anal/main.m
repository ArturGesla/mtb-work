mu=0.04;;
r=sqrt(mu);
t=0;
ta=[];
eva=[];
%%
for i=1:32
t=t+0.1;
ta=[ta;t];
% r=0;
g=1; %gamma
u=r*cos(t);
v=r*sin(t)/g;
w=r^2;

f=[mu*u-v*g-w*u;
    mu*v+u/g-v*w;
    -w+u^2+v^2*g*g];
J=[mu-w, -g,-u;
    1/g,mu-w,-v;
    2*u,2*v*g*g,-1];
eig(J)
P=eye(3)-f*f'/norm(f)/norm(f)
Jp=P*J*P
[evc,ev]=eig(Jp)
eva=[eva,diag(ev)];
end
%%
plot(mod(ta,0),eva','-x')
