clf;
a1=load("evsCheb.mat");

% loglog(a1.narr,abs(a1.uarr-a1.uarr(end,:)),'x-'); hold on;

a2=load("evsFD.mat");
loglog(a2.narr,abs(a2.uarr-a1.uarr(end,:)),'x-'); hold on;

%%

% cos(pi/2*xc)*cos(pi/2*yc)
clf; leg=[];
set(gcf,"Position",   [515   164   807   618])

a=load("centerFD2.mat"); 
loglog(a.narr,abs(a.uarr-2/pi/pi),'x-'); hold on; leg=[leg;"finite diff 2nd order"]; 
loglog(a.narr,a.narr.^(-2),'-'); leg=[leg;"n^{-2} slope"]; 
loglog(a.narr,a.narr.^(-4),'-'); leg=[leg;"n^{-4} slope"]; 
a=load("centerFD4.mat"); 
loglog(a.narr,abs(a.uarr-2/pi/pi),'x-'); hold on; leg=[leg;"finite diff compact 4th order"]; 

a=load("chebcenter.mat"); 
loglog(a.narr,abs(a.uarr-2/pi/pi),'x-'); hold on; leg=[leg;"Chebyshev"]; 
xlabel("n,number of points in one direction"); ylabel("error=numerical solution at (0,0)-2/pi/pi");
title("numercial solution of Poisson with f=cos(pi/2*xc)*cos(pi/2*yc) in [-1,1]x[-1,1] square");
legend(leg); grid on;