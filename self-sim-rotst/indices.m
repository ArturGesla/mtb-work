iip=(iz-1)*4+1;
iipzp=iip+4;
iif=(iz-1)*4+2;
iifzp=iif+4;
iifzm=iif-4;
iig=(iz-1)*4+3;
iigzp=iig+4;
iigzm=iig-4;
iih=(iz-1)*4+4;
iihzp=iih+4;
iihzm=iih-4;

iik=Nz*4+1;

if (iz>1 && iz<Nz) 
hz=(zw(iz)-zw(iz-1)); 
hzw=(zc(iz+1)-zc(iz));
    
fn=((zw(iz-1)-zc(iz-1))*u(iif)+(zc(iz)-zw(iz-1))*u(iifzm))/(zc(iz)-zc(iz-1));
fs=((zw(iz)-zc(iz))*u(iifzp)+(zc(iz+1)-zw(iz))*u(iif))/(zc(iz+1)-zc(iz));

azum=1/(zc(iz)-zc(iz-1))/hz;
azup=1/(zc(iz+1)-zc(iz))/hz;

gn=((zw(iz-1)-zc(iz-1))*u(iig)+(zc(iz)-zw(iz-1))*u(iigzm))/(zc(iz)-zc(iz-1));
gs=((zw(iz)-zc(iz))*u(iigzp)+(zc(iz+1)-zw(iz))*u(iig))/(zc(iz+1)-zc(iz));

azvm=1/(zc(iz)-zc(iz-1))/hz;
azvp=1/(zc(iz+1)-zc(iz))/hz;

hn=((zc(iz)-zw(iz-1))*u(iih)+(zw(iz)-zc(iz))*u(iihzm))/(zw(iz)-zw(iz-1));
hs=((zc(iz+1)-zw(iz))*u(iihzp)+(zw(iz+1)-zc(iz+1))*u(iih))/(zw(iz+1)-zw(iz));

azwm=1/(zw(iz)-zw(iz-1))/hzw;
azwp=1/(zw(iz+1)-zw(iz))/hzw;

pn=u(iip);
ps=u(iipzp);

end