clear;
om1=1;
om2=0;
om3=1;

s1=om1/om1; s2=om2/om1; s3=om3/om1;

[R,Z]=meshgrid(0:0.01:1,0:0.02:1);
h=0.1; l=1; delta=h/l;
lamN=@(n) n*pi/delta;

%%
v1=((s2-s1)*Z+s1).*R;
v2=v1*0;
N=200;
dvarr=[];
for n=1:N
    arg=lamN(n)*R;


    % An=@(n) 2/n/pi/besseli(1,lamN(n),0)*((s2-s3)*-1^n+(s3-s1));
    %     dv=An(n)*besseli(1,arg,0).*sin(n*pi*Z);

    An=@(n) 2./n./pi./besseli(1,lamN(n),1).*((s2-s3)*((-1).^n)+(s3-s1));
    dv=An(n)*besseli(1,arg,1).*sin(n*pi*Z).*exp(arg-lamN(n));
dvarr=[dvarr,dv(51,end)];

    v2=v2+dv;
end
v=v1+v2;
%
mesh(R,Z,v)
% pcolor(R,Z,v1); colorbar(); shading interp; colormap(parula(8))

%%

n=1:100;
% semilogy(n,abs(An(n).*besseli(1,lamN(n),1).*sin(n*pi*0.5)))
loglog(n,abs(An(n).*besseli(1,lamN(n),1)));
hold on;
loglog(n,1./n)

%% ss reg
vss=1-Z(:,2:end);
vst=v(:,2:end)./R(:,2:end);
sserr=vecnorm(vst-vss)./vecnorm(vss);
semilogy(R(1,2:end),sserr)

%%
clf; set(gcf,"Position",[         256         333        1009         389])
subplot(3,2,1);
a=load("shroudrot.mat");
pcolor(a.R*10,a.Z,a.v*10); colorbar(); shading interp; colormap(parula(10))
caxis([0 10]); title("$u_{\theta}$ shroud rotating, $Re\rightarrow 0$");
fnts=10; jfm_plt_aid_comm;

subplot(3,2,3);
a=load("shroudfix.mat");
pcolor(a.R*10,a.Z,a.v*10); colorbar(); shading interp; colormap(parula(10))
caxis([0 10]); title("$u_{\theta}$ shroud fixed, $Re\rightarrow 0$");
fnts=10; jfm_plt_aid_comm;

subplot(3,2,5);
a=load("shroudrot.mat");
pcolor(a.R*10,a.Z,a.v1*10); colorbar(); shading interp; colormap(parula(10))
caxis([0 10]); title("$u_{\theta}$ shroud lin, $Re\rightarrow 0$");
fnts=10; jfm_plt_aid_comm;

subplot(3,2,[2,4,6]);
a=load("shroudrot.mat");
vss=1-Z(:,2:end);
vst=v(:,2:end)./R(:,2:end);
sserr=vecnorm(vst-vss)./vecnorm(vss);
semilogy(R(1,2:end),sserr,'-o')

hold on; grid on;
a=load("shroudrot.mat");
vss=1-Z(:,2:end);
vst=v1(:,2:end)./R(:,2:end);
sserr=vecnorm(vst-vss)./vecnorm(vss);
semilogy(R(1,2:end),sserr)


a=load("shroudfix.mat");
vss=1-Z(:,2:end);
vst=v(:,2:end)./R(:,2:end);
sserr=vecnorm(vst-vss)./vecnorm(vss);
semilogy(R(1,2:end),sserr,'-x')

title("error of non self similarity, $Re\rightarrow 0$"); 
legend("shround rot","shround lin","shroud fix","Location","best")
text(0.1,1e-3,"error = 0.01 at around r=8.55")
fnts=12; jfm_plt_aid_comm;
