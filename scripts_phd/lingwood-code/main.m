%% base flow
n=20;
neq=3;
u=zeros(3*n,1);
zmax=20;
dz=zmax/n;
z=0:dz:zmax-dz;

triplets=[];

g=zeros(3*n,1);

%parms
Ro=-1;
Co=2-Ro-Ro^2;

%BC
fBottom=0;
gBottom=0;
hBottom=0;

fTop=0;
gTop=1;

%bottom
i=1;
iif=i;
iig=iif+1;
iih=iig+1;

g(iif)=u(iif)-fBottom;
g(iig)=u(iig)-gBottom;
g(iih)=u(iih)-hBottom;

%main body
for i=2:n-1
iif=i;
iig=iif+1;
iih=iig+1;
iifzp=iif+neq;
iifzm=iif-neq;

g(iif)=Ro*(u(iif)^2+u(iih)*(u(iifzp)-u(iifzm))/2/dz-(u(iig)^2-1))-Co*(u(iig)-1)-(u(iifzp)-2*u(iif)+u(iifzm))/dz/dz;
g(iig)=u(iig)-gBottom;
g(iih)=u(iih)-hBottom;

end
