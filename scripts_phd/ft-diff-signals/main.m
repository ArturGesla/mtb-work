clf
%%
x=0:0.01:1; x=x(1:end-1);
% y=ones(length(x),1); y(51:end)=0;
% y=x*2; y(51:end)=-2*x(51:end)+2;
y=ones(length(x),1); y(51:end)=0; 

%
set(gcf,"Position",[76         369        1436         413])
subplot(1,2,1); hold on;
plot(x,y); grid on;

subplot(1,2,2); 
z=abs(fft(y));
loglog(z(1:end/2),'x'); grid on; ylim([1e-3 1e2]); hold on;
