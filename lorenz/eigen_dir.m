sigma=10;
b=8/3;
r=25;

c=sqrt(b*(r-1));
x=c;
y=c;
z=r-1;

x=0;
y=0;
z=0;

J=[-sigma, sigma,0;r-z,-1,-x;y,x,-b];

[v,u]=eig(J);
u=diag(u);

u
v
%%
U=[];
R=[];

for r=0.1:0.1:4
J=[-sigma, sigma,0;r-z,-1,-x;y,x,-b];

[v,u]=eig(J);
u=diag(u);

U=[U,u];
R=[R,r];
end

plot(R,U')