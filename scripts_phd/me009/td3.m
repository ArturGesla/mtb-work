nx=5;
ny=4;
m=diag(ones(nx,1))*4-diag(ones(nx-1,1),-1)-diag(ones(nx-1,1),1);
n=eye(nx)*-1;
A=[];

i=1;
A=[A;repmat(zeros(nx),1,i-1-2), m,n,repmat(zeros(nx),1,ny-i-2+1)];
for i=2:ny-1
    A=[A;repmat(zeros(nx),1,i-2), n,m,n,repmat(zeros(nx),1,ny-i-2+1)];
end
i=ny;
A=[A;repmat(zeros(nx),1,i-2), n,m,repmat(zeros(nx),1,ny-i-2+1)];

gs=0;
gw=1;
gn=0; 
ge=0;

b=[ones(nx,1);repmat(zeros(nx,1),ny-1,1)]*gs+repmat([1;zeros(nx-1,1)],ny,1)*gw+repmat([zeros(nx-1,1);1],ny,1)*ge+gn*[repmat(zeros(nx,1),ny-1,1);ones(nx,1);];

%% solve

x=A\b;
xx=reshape(x,nx,ny)';
 contour(xx);

 %%
 evs=eig(A);
 plot(real(evs),imag(evs),'x')
