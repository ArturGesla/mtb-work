clc; clear;
% N=3;
% beta=0.5./sqrt(1-(2*(1:N)).^(-2))
% T=diag(beta,1)+diag(beta,-1) % rec relation coeffs
% [V,D]=eig(T)
% x=diag(D); [x,i]=sort(x);
% w=2*V(1,i).^2;

%

varr=[];
varr2=[];

f=@(x) cos(x*pi/2);  vext=4/pi; % not periodic 
% f=@(x) 1+cos(x*pi);  vext=; % int of constant wtf 

for i=1:25
N=i;
[w,x]=quadraturewx(N);
v=f(x)'*w';
varr=[varr;v];
xt=-1:2/N:1;
varr2=[varr2;trapz(xt,f(xt))];
end

%%
clf;
loglog(abs(varr-4/pi)); hold on
loglog(abs(varr2-4/pi)); hold on
loglog(5*[1:length(varr2)].^(-2))