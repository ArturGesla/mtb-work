
cd /home/gesla/Documents/git/rotst2/build;
N=3; name="";
a=importdata("u-"+num2str(N)+".dat"+name);
x=reshape(a(:,1),[N+1,N]);
y=reshape(a(:,2),[N+1,N]);
u=reshape(a(:,3),[N+1,N]);
mesh(x,y,u);

xp=x/2+1/2;
yp=y/2+1/2;
% mesh(x,y,u -3*(xp.^3.*yp.^2));

%%

a=importdata("v-"+num2str(N)+".dat"+name);
x=reshape(a(:,1),[N+1,N+1]);
y=reshape(a(:,2),[N+1,N+1]);
v=reshape(a(:,3),[N+1,N+1]);
mesh(x,y,v);
title("v"); 

%%

a=importdata("w-"+num2str(N)+".dat"+name);
x=reshape(a(:,1),[N,N+1]);
y=reshape(a(:,2),[N,N+1]);
w=reshape(a(:,3),[N,N+1]);
% mesh(x,y,w);
xp=x/2+1/2;
yp=y/2+1/2;
mesh(x,y,w +4*(xp.^2.*yp.^3));
title("w"); 


%%

a=importdata("p-"+num2str(N)+".dat"+name);
x=reshape(a(:,1),[N-1,N-1]);
y=reshape(a(:,2),[N-1,N-1]);
p=reshape(a(:,3),[N-1,N-1]);
mesh(x,y,p);
title("p"); 

%%

% (chebInt(1)+chebIntX(1))*(chebInt(1)) %1
1/8*(chebInt(1)+chebIntX(1))*(1/2*(chebInt(2)+chebInt(0)))   %2, 1/18, ok