function dydt=subhopf(t,y,lam)
dydt=zeros(1,1);
dydt=[y(1)*((lam-1)/lam-y(1))];
end