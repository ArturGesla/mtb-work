clc; clear; close all;



%%
xa=[x];
for i=1:10
    f=xa(end)^2-3;
    jac=2*xa(end);
    xa(end+1)=xa(end)-f/jac;
end

%%
x=20;
xa=[x];
for i=1:3000
    f=xa(end)^2-3;
    jacT=2*xa(end);
    xa(end+1)=xa(end)-jacT*f*0.001;
    xa(end)
end

%%
x=-2.5:0.1:2.5;
y=x.*(x-2).*(x+2)+3.5;
y=x.^3-4.*x+3.5;
plot(x,y);
grid on;
%%
hold on; plot(xa(end),xa(end).^3-4.*xa(end)+3.5,'x');
%%
x=2; 
xa=[x];
for i=1:10000
    f=xa(end).^3-4.*xa(end)+3.5;
    jac=3*xa(end)^2-4;
    % xa(end+1)=xa(end)-f/jac;
    xa(end+1)=xa(end)-jac*f*1e-2;
    xa(end)
end
plot(xa)
