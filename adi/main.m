lx=3;
ly=2;
nx=25;
ny=20;
dt=0.01;

hx=lx/(nx-2);
hy=ly/(ny-2);
u=ones(nx*ny,1)*0-1;
u=reshape(u,[ny,nx]);
u(2:end-1,2:end-1)=1;
u=reshape(u,[ny*nx,1]);
%%
L=zeros(nx*ny);
B=L;
rhs=zeros(nx*ny,1);
for i=2:nx-1
    for j=2:ny-1
        ii=(i-1)*ny+j;
        L(ii,ii)=-2/hx/hx-2/hy/hy-1/dt;
        L(ii,ii+ny)=1/hx/hx;
        L(ii,ii-ny)=1/hx/hx;
        L(ii,ii+1)=1/hy/hy;
        L(ii,ii-1)=1/hy/hy;
        rhs(ii)=(u(ii-1)-2*u(ii)+u(ii+1))/hy/hy+(u(ii-ny)-2*u(ii)+u(ii+ny))/hx/hx;
        B(ii,ii)=1;
    end
end
%corners
i=1;j=1; ii=(i-1)*ny+j; L(ii,ii)=1; rhs(ii)=0;
i=nx;j=1; ii=(i-1)*ny+j; L(ii,ii)=1; rhs(ii)=0;
i=1;j=ny; ii=(i-1)*ny+j; L(ii,ii)=1; rhs(ii)=0;
i=nx;j=ny; ii=(i-1)*ny+j; L(ii,ii)=1; rhs(ii)=0;
%BC
for i=1
    for j=2:ny-1
        ii=(i-1)*ny+j;
        L(ii,ii)=1;
        L(ii,ii+ny)=1;
        rhs(ii)=0;
    end
end
for i=nx
    for j=2:ny-1
        ii=(i-1)*ny+j;
        L(ii,ii)=1;
        L(ii,ii-ny)=1;
        rhs(ii)=0;
    end
end
for i=2:nx-1
    for j=1
        ii=(i-1)*ny+j;
        L(ii,ii)=1;
        L(ii,ii+1)=1;
        rhs(ii)=0;
    end
end
for i=2:nx-1
    for j=ny
        ii=(i-1)*ny+j;
        L(ii,ii)=1;
        L(ii,ii-1)=1;
        rhs(ii)=0;
    end
end

L;
%  plot(eig(L),'x')
%

% OP=diag(1/dt*ones(nx*ny,1))-L;
% du=OP\(rhs);
du=L\(-rhs);
u=u+du;
%
up=reshape(u,[ny,nx]);
mesh(up);
zlim([0,1]);
%%
[v,ev]=eig(L,B);
%%
iev=iev+1;
up=reshape(v(:,iev),[ny,nx]);
ev(iev,iev)
mesh(up);