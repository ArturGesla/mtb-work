clear;
om1=1;
om2=0;
om3=1;

s1=om1/om1; s2=om2/om1; s3=om3/om1;

[R,Z]=meshgrid(0:0.01:1,0:0.01:1);
h=1; l=1; delta=h/l;
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
mesh(R,Z,v2)
% pcolor(R,Z,v2); colorbar(); shading interp; colormap(parula(8))

%%

n=1:100;
% semilogy(n,abs(An(n).*besseli(1,lamN(n),1).*sin(n*pi*0.5)))
loglog(n,abs(An(n).*besseli(1,lamN(n),1)));
hold on;
loglog(n,1./n)