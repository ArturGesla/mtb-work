function dydt=subhopfp_r3(t,y,mu,b,om,vG)
dydt=zeros(2,1);
dydt=[y(1)*(mu+y(1)^2+y(2)^2)-y(2)*(om+b*(y(1)^2+y(2)^2)); y(2)*(mu+y(1)^2+y(2)^2)+y(1)*(om+b*(y(1)^2+y(2)^2))]+(vG);
end