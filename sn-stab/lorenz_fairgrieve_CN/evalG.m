function g=evalG(x,y,z,sigma,r,b)
g=zeros(3,1);
g(1)=sigma*(y-x);
g(2)=r*x-y-x*z;
g(3)=x*y-b*z;

end
% gstab(ix)=(sigma*(u(iy)-u(ix)))      ;
%     gstab(iy)=(r*u(ix)-u(iy)-u(ix)*u(iz))   ;
%     gstab(iz)=(u(ix)*u(iy)-b*u(iz))           ;