function dydt=subhopfw(t,y,mu,b,om,lamc,k,ytau)
dydt=zeros(3,1);
dydt=[mu*y(1)+y(1)^3-y(1)^5+y(3)*y(1); om+b*y(1)^2;lamc*y(3)-k*(y(1)-ytau(1))*y(1)];
end