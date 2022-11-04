r=28;
lorenz = @(t,x) [10*(x(2)-x(1));r*x(1)-1*x(1)*x(3)-x(2);x(1)*x(2)-(8/3)*x(3)]; % Anonymous Function
[T,X] = ode45(lorenz, [0 :1.5586/(np):1.5586], [-6.2262, -11.0027,13.0515]);
%%
% hold on;
%  plot3(X(:,1),X(:,2),X(:,3))
%  plot3(X(1,1),X(1,2),X(1,3),'ro')
%  disp(sqrt(8/3*(r-1)))
%  %%
% plot(mod(T,1.5586),X,'x')