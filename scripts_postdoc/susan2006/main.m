a=importdata("data2");
c=a(:,4:11);
%%
clf;
% i=5;
for i=1:17
om0=c(i,1);
om1=c(i,2);
om2=c(i,3);
u0=c(i,4);
u1=c(i,5);
u2=c(i,6);
r1=c(i,7);
r2=c(i,8);

%

r=-1:0.0021:1;

uz=u0+u1*exp(-r.^2/r1^2)+u2*exp(-r.^2/r2^2);
uth=om0*r+om1*r1^2./r.*(1-exp(-r.^2/r1^2))+om2*r2^2./r.*(1-exp(-r.^2/r2^2));


up=uth;
if i==3
plot(r,up,'r',LineWidth=3)
else
plot(r,up)
end

hold on; grid on;
end