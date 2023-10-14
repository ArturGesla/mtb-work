function dydt=langford(t,y,lam,c,k,vG)
dydt=zeros(3,1);
l=lam;
dydt=[(l-3)*y(1)-1/4*y(2)+y(1)*(y(3)+0.2*(1-y(3)^2));
    (l-3)*y(2)+1/4*y(1)+y(2)*(y(3)+0.2*(1-y(3)^2));
    l*y(3)-(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))];
end