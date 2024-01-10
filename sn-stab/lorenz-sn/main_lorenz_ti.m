% r=28; %r=320; 
nP=1;
lorenz = @(t,x) [10*(x(2)-x(1));r*x(1)-1*x(1)*x(3)-x(2);x(1)*x(2)-(8/3)*x(3)]; % Anonymous Function
% [T,X] = ode45(lorenz, [0 :T/(np):1*T], [-6.2262, -11.0027,13.0515]);
T=0.6779; T= 0.6803; [t,X] = ode45(lorenz, [0 :T/(np-2):T/(np-2)*(np-1)*nP], [ 9.2898   11.1031   22.2085]);
% T=1/0.868; [t,X] = ode45(lorenz, [0 :T/(np):1*T], [-5.6775    2.1838  134.4854]);
% plot3(X(:,1),X(:,2),X(:,3))
plot(X(:,1),X(:,2))
%%
% hold on;
% last=1000;
%  plot3(X(:,1),X(:,2),X(:,3))
%  plot3(X(1,1),X(1,2),X(1,3),'ro')
%  disp(sqrt(8/3*(r-1)))
% %  %%
% % plot(mod(T,1.5586),X,'x')
% 
% %% fft
% x=X(:,1);
% df=1/(t(end)-t(1));
% dt=mean(diff(t));
% f=0:df:1/dt;
% z=abs(fft(x));
% plot(f,z)
% %%
% ff=0.868;
% T=1/ff;
