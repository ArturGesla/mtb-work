lorenz = @(t,x) [10*(x(2)-x(1));r*x(1)-1*x(1)*x(3)-x(2);x(1)*x(2)-(8/3)*x(3)]; % Anonymous Function
% [T,X] = ode45(lorenz, [0 :T/(np):1*T], [-6.2262, -11.0027,13.0515]);
T=0.6779; T= 0.6803; [t,X] = ode45(lorenz, [0 :T/(np):150*T], [ 9.2898   11.1031   22.2085]);
% T=1/0.868; [t,X] = ode45(lorenz, [0 :T/(np):1*T], [-5.6775    2.1838  134.4854]);
% plot3(X(:,1),X(:,2),X(:,3))
plot(X(:,1),X(:,2))
%%
close all;

en=sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2);
plot(t,en);

%% not so good
