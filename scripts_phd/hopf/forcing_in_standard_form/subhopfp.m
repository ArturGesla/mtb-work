function dydt=subhopfp(t,y,mu,b,om,G)
dydt=zeros(2,1);
dydt=[mu*y(1)+y(1)^3-y(1)^5; om+b*y(1)^2]+diag(G);
end