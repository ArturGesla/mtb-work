addpath     '/home/gesla/Documents/git/rotst2/scripts/source_for_mtb'
cd rotst ;
% a=importdata("u10-64-74.dat");
% a=loadrs2(a,64,74);
% a=importdata("u10-151-151.dat");
% a=loadrs2(a,151,151);
a=importdata("u10-301-301.dat");
a=loadrs2(a,301,301);
cd ..;
% mesh(a.xc,a.zc,a.u);
% mesh(a.xc,a.zc,a.w);
%
% u interpl

uC=interp2(a.xu(1,1:end)*2/R-1,a.zc(:,1)*2/H-1,a.u(1:end,1:end-1),ru,zu');
wC=interp2(a.xc(1,1:end)*2/R-1,a.zw(:,1)*2/H-1,a.w(1:end-1,1:end),rw,zw');

% mesh(ru,zu',uC);
% mesh(rw,zw',wC);
%
%
% tic;
% A=zeros(length(ru)*length(zu));
A=zeros((N+1)^2);
for ix=1:length(ru)
    for iy=1:length(zu)
        ip=iy+(ix-1)*(N+1); xc=ru(ix); zc=zu(iy);
        for ikx=0:length(ru)-1
            for iky=0:length(zu)-1
                ikp=iky+1+(ikx-1+1)*(N+1);
                A(ip,ikp)=A(ip,ikp)+Tn(xc,ikx)*Tn(zc,iky);
            end
        end
    end
end
for ix=N+1
    for iy=1:length(zu)
        ip=iy+(ix-1)*(N+1);
        A(ip,ip)=1;

    end
end
%
au=A\[reshape(uC,[length(ru)*length(zu),1]);zeros(N+1,1)];
%

% mesh(reshape(log(abs(au)),[N+1,N+1]))
% mesh(reshape(log(abs(aw)),[N+1,N+1]))
%%
A=zeros((N+1)^2);
for ix=1:length(rw)
    for iy=1:length(zw)
        ip=iy+(ix-1)*(N+1); xc=rw(ix); zc=zw(iy);
        for ikx=0:length(rw)-1
            for iky=0:length(zw)-1
                ikp=iky+1+(ikx-1+1)*(N+1);
                A(ip,ikp)=A(ip,ikp)+Tn(xc,ikx)*Tn(zc,iky);
            end
        end
    end
end
for ix=1:length(rw)
    for iy=N+1
        ip=iy+(ix-1)*(N+1); 
         A(ip,ip)=1;
    end
end
%
aw=A\[reshape([wC;zeros(1,N+1)],[(N+1)^2,1]);];
% toc;