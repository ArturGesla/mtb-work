addpath     '/home/gesla/Documents/git/rotst2/scripts/source_for_mtb'
cd rotst ;
% a=importdata("u10-64-74.dat");
% a=loadrs2(a,64,74);
% a=importdata("u10-151-151.dat");
% a=loadrs2(a,151,151);
a=importdata("u10-301-301.dat");
% a=importdata("u1-301-301.dat"); %smoother
% a=importdata("u500-151-151.dat"); %smoother
a=loadrs2(a,301,301);
% a=loadrs2(a,151,151);
cd ..;
% mesh(a.xc,a.zc,a.u);
% mesh(a.xc,a.zc,a.w);
%
% u interpl

pC=interp2(a.xc(1,1:end)*2/R-1,a.zc(:,1)*2/H-1,a.p(1:end,1:end),rp,zp');
uC=interp2(a.xu(1,1:end)*2/R-1,a.zc(:,1)*2/H-1,a.u(1:end,1:end-1),ru,zu');
vC=interp2(a.xc(1,1:end)*2/R-1,a.zc(:,1)*2/H-1,a.v(1:end,1:end),rv,zv'); vC(end,end)=1; 
% vC(1,end)=1;
wC=interp2(a.xc(1,1:end)*2/R-1,a.zw(:,1)*2/H-1,a.w(1:end-1,1:end),rw,zw');

%  mesh(ru,zu',uC);
% mesh(rw,zw',wC);
% mesh(rv,zv',vC);
mesh(rp,zp',pC);
%%
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
%% v map

A=zeros((N+1)^2);
for ix=1:length(rv)
    for iy=1:length(zv)
        ip=iy+(ix-1)*(N+1); xc=rv(ix); zc=zv(iy);
        for ikx=0:length(rv)-1
            for iky=0:length(zv)-1
                ikp=iky+1+(ikx-1+1)*(N+1);
                A(ip,ikp)=A(ip,ikp)+Tn(xc,ikx)*Tn(zc,iky);
            end
        end
    end
end
%
% for ix=1:length(rw)
%     for iy=N+1
%         ip=iy+(ix-1)*(N+1); 
%          A(ip,ip)=1;
%     end
% end
%
av=A\[reshape([vC],[(N+1)^2,1]);];
%% v remap
vv=zeros(length(zv),length(rv));
%
for ix=1:length(rv)
    for iy=1:length(zv)
        ip=iy+(ix-1)*(N+1); xc=rv(ix); zc=zv(iy);
        for ikx=0:length(rv)-1
            for iky=0:length(zv)-1
                ikp=iky+1+(ikx-1+1)*(N+1);
                vv(iy,ix)=vv(iy,ix)+av(ikp)*Tn(xc,ikx)*Tn(zc,iky);
            end
        end
    end
end
%
mesh(zv,rv,vv)
%% p map
A=zeros((N+1)^2);
%
for ix=1:length(rp)
    for iy=1:length(zp)
        ip=(iy+1)+(ix-1+1)*(N+1); xc=rp(ix); zc=zp(iy);
        for ikx=0:length(rp)-1
            for iky=0:length(zp)-1
                ikp=iky+1+1+(ikx-1+1+1)*(N+1);
                A(ip,ikp)=A(ip,ikp)+Tn(xc,ikx)*Tn(zc,iky);
            end
        end
    end
end
%
for ix=1:length(rv)
    for iy=1:length(zv)
        ip=iy+(ix-1)*(N+1);
        if (ix==1)
         A(ip,ip)=1;
        end
        if (ix==length(rv))
         A(ip,ip)=1;
        end
        if (iy==1)
         A(ip,ip)=1;
        end
        if (iy==length(zv))
         A(ip,ip)=1;
        end

    end
end
%
ppC=[zeros(1,N+1);[zeros(N-1,1),pC,zeros(N-1,1)];zeros(1,N+1)];
%
ap=A\[reshape([ppC],[(N+1)^2,1]);];
