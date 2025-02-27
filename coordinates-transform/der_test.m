clc; clear;
f=@(r) r^2;

eps=1e-2;
df=@(r) (f(r+eps)-f(r-eps))/2/eps;
d2f=@(r) (f(r+eps)+f(r-eps)-2*f(r))/eps/eps;

z=1; L=2;
% fth=@(th) ((L-z)*tan(th)).^2;
fth=@(th,z) ((L-z)*tan(th));
% dfth=@(th) (fth(th+eps)-fth(th-eps))/2/eps*cos(th)^2/(L-z); %original
dfth=@(th,z) (fth(th,z+eps)-fth(th,z-eps))/2/eps+(fth(th+eps,z)-fth(th-eps,z))/2/eps*cos(th)*sin(th)/(L-z);
% d2fth=@(th) -2*sin(th)*cos(th)/(L-z)*(fth(th+eps)-fth(th-eps))/2/eps*cos(th)^2/(L-z)+(cos(th)^2/(L-z))^2*(fth(th+eps)+fth(th-eps)-2*fth(th))/eps/eps;
% d2fth=@(th) cos(th)^2/(L-z)*(-2*sin(th)*cos(th)/(L-z)*(fth(th+eps)-fth(th-eps))/2/eps+cos(th)^2/(L-z)*(fth(th+eps)+fth(th-eps)-2*fth(th))/eps/eps);
d2fthdy=@(th,z) (fth(th,z+eps)+fth(th,z-eps)-2*fth(th,z))/eps/eps+sin(th)*cos(th)/(L-z)^2*(fth(th+eps,z)-fth(th-eps,z))/2/eps+sin(th)*cos(th)/(L-z)*((fth(th+eps,z+eps)-fth(th-eps,z+eps))/2/eps-(fth(th+eps,z-eps)-fth(th-eps,z-eps))/2/eps)/2/eps;
% +cos(th)^2/(L-z)*(-2*sin(th)*cos(th)/(L-z)*(fth(th+eps)-fth(th-eps))/2/eps+cos(th)^2/(L-z)*(fth(th+eps)+fth(th-eps)-2*fth(th))/eps/eps);

%%
% d2fth=

tha=[0:0.01:1];
va=tha*0;
for i=1:length(tha)
% %     va(i)=d2fth(tha(i));
    va(i)=d2fthdy(tha(i),0.1);
%     va(i)=dfth(tha(i),0);
end
%
plot(tha,va)