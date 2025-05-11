function dydt=subhopfp_r3(t,y,mu,b,om,vG)
dydt=zeros(3,1);
dydt=[(y(1)*mu-y(2)*b-y(1)*y(3)); (y(2)*mu+y(1)/b-y(2)*y(3)); -y(3)+y(1)*y(1)+y(2)*y(2)*b^2]+(vG);
end