function dydt=subhopfp_r3(t,y,m,c,k,vG)
dydt=zeros(3,1);
dydt=[y(2); 1/m*(-c*y(2)-k*y(1)); 0]+(vG);
end