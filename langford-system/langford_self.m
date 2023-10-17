
l=1.8; close all;
langford = @(t,y) [(l-3)*y(1)-1/4*y(2)+y(1)*(y(3)+0.2*(1-y(3)^2));
    (l-3)*y(2)+1/4*y(1)+y(2)*(y(3)+0.2*(1-y(3)^2));
    l*y(3)-(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))]; % Anonymous Function
 [t,X] = ode45(langford, [0,1000], [ 0.2334 -0.9401 1.0858]);
% plot(X(:,1),X(:,2))
plot3(X(:,1),X(:,2),X(:,3))

%%
close all;
figure("Position",[1921           1        1920         997])
for i=1:10
subplot(2,5,i)
% l=2.01; 
% l=1.98+i/10*(2.023-1.98); 
l=1.5+i/10*(1.98-1.5); 
langford = @(t,y) [(l-3)*y(1)-1/4*y(2)+y(1)*(y(3)+0.2*(1-y(3)^2));
    (l-3)*y(2)+1/4*y(1)+y(2)*(y(3)+0.2*(1-y(3)^2));
    l*y(3)-(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))]; % Anonymous Function
 [t,X] = ode45(langford, [0,1000], [ 0.2334 -0.9401 1.0858]);
% plot(X(:,1),X(:,2))
plot3(X(:,1),X(:,2),X(:,3))
xlim([-1.2 1.2]); ylim([-1.2 1.2]); zlim([0 2]);
title("lam: "+num2str(l));
grid on;
end
%%
exportgraphics(gcf,"plot1.png","Resolution",300);

%%
x=X(end-1000:end,1);
y=X(end-1000:end,2);
z=X(end-1000:end,3);
t1=t(end-1000:end);
% plot3(x,y,z)
plot(t1,x); clf;
[f,z]=ft(t1,x); semilogy(f,z); hold on; grid on;
[f,z]=ft(t1,y); semilogy(f,z); hold on; grid on;
[f,z]=ft(t1,z); semilogy(f,z); hold on; grid on;
xlim([0 0.2]);
ylim([1e-6 1])