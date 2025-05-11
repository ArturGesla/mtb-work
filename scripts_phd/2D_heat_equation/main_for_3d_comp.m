clc; clear;
uarr=[];
%
nx=256*4; k=1;
narr=[];

x=linspace(-1,1,nx)'; dx=mean(diff(x));x=-dx/2:dx:1+dx/2; nx=length(x);

u=rand(nx*nx,1); 
f=rand(nx*nx,1);
% f=f-mean(f); 
% f=ones(nx*nx,1); 
[g,jac]=evalJacAndRhs(x,u,f,k);

u=u-jac\g;

[g,jac]=evalJacAndRhs(x,u,f,k);
norm(g)
[U,S,V]=svds(jac,5,"smallest");
%%

up=u;
% up=V(:,end-2);
mesh(x,x,reshape(up,[nx,nx]))
% h=pcolor(x,x,reshape(u,[nx,nx])); 
% shading interp; colormap(parula(16))
% h.EdgeAlpha=0; colormap(parula(16))
% mesh(x,x,reshape(g,[nx,nx]))

%%
% h=pcolor(x,x,reshape(f,[nx,nx])); 
h=pcolor(x,x,reshape(up,[nx,nx])); 
% shading interp; colormap(parula(16))
h.EdgeAlpha=0; colormap(parula(16))
