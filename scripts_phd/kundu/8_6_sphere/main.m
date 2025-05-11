clc; clear; 
x=-2:0.01:2;
y=x;

psi=zeros(length(y),length(x));
U=1;
a=1;

for ix=1:length(x)
    for iy=1:length(y)
        xc=x(ix); yc=y(iy); th=atan(yc/xc); r=sqrt(xc^2+yc^2);
       
        psi(iy,ix)=U*r^2*sin(th)^2*(1/2-3*a/4/r+a^3/4/r^3);
%         psi(iy,ix)=th;
    end
end

%
[c,h]=contour(x,y,psi,[0:0.1:1].^4); colorbar(); hold on; %clabel(c,h); 
% plot(x,real(sqrt(a^2-y.^2)),':k');
% plot(x,-real(sqrt(a^2-y.^2)),':k');
colormap(hot(10));