clc; clear;
r=24; sigma=10; beta=8/3;





x0=[5;5;5];
rng(1);
x0=[sqrt(beta*(r-1));-sqrt(beta*(r-1));r-1]+rand(3,1)*0.1;
n=2000;
X=zeros(3,n);
X(:,1)=x0;



fa=zeros(n-1,1);

for i=1:n-1
    x=X(1,i); y=X(2,i); z=X(3,i); 
f=[sigma*(y-x); x*(r-z)-y; x*y-beta*z];
jac=[-sigma, sigma, 0; r-z, -1, -x; y,x,-beta];
X(:,i+1)=X(:,i)-jac\f;
% X(:,i+1)=X(:,i)-jac'*f*1e-3;
fprintf("it\t%d\tnorm\t%e\n",i,norm(f));
fa(i)=norm(f);
if(norm(f)<1e-10)
    break;
end
end


%%


% xc=-10:0.1:10; yc=-9:0.1:9; zc=zeros(length(xc),length(yc));
xc=-15:0.1:15; yc=-14:0.1:14; yc=[0]; zc=zeros(length(xc),length(yc));
% xc=-15:1:15; yc=-14:0.1:14; zc=zeros(length(xc),length(yc));
% xc=-150:1:150; yc=-140:0.1:140; zc=zeros(length(xc),length(yc));
% dx=0.1; yc=1:dx:5; xc=-10:dx:-6; zc=zeros(length(xc),length(yc));

for ix=1:length(xc)
    for iy=1:length(yc)

x0=[xc(ix);yc(iy);10];
% rng(1);
% x0=[sqrt(beta*(r-1));-sqrt(beta*(r-1));r-1]*0+rand(3,1)*0.1;
n=20000;
X=zeros(3,n);
X(:,1)=x0;





for i=1:n-1
    x=X(1,i); y=X(2,i); z=X(3,i); 
f=[sigma*(y-x); x*(r-z)-y; x*y-beta*z];
jac=[-sigma, sigma, 0; r-z, -1, -x; y,x,-beta];
% X(:,i+1)=X(:,i)-jac\f;
X(:,i+1)=X(:,i)-jac'*f*1e-3;

% fprintf("it\t%d\tnorm\t%e\n",i,norm(f));
if(norm(f)<1e-3)
    break;
else
    error("no conv");
end
end
zc(ix,iy)=X(1,i);

    end
    ix
end

%%
pcolor(yc,xc,zc); shading interp; axis equal; colorbar(); colormap(parula(3)); clim([-8 8])