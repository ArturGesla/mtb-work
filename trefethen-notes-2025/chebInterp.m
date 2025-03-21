n=41;
x=-cos((0:n)/n*pi);
% x=-1:2/(n):1;
% y1=x;
% y1=@(x) exp(x);
% y1=@(x) sin(pi*(x-0.1));
% y1=@(x) 1./(1+25*x.^2);
% y1=@(x) 1./(1+x.^2);
% y1=@(x) abs(x);
y1=@(x) sin(x*pi);
a=[];
%% cheb basis
T=cos([0:n]'*acos(x));
% a=T\exp(x')
a=T'\(y1(x)'); cond(T')
x2=-1:2/n/2:1;
y2=a'*cos([0:n]'*acos(x2));
a(end)
%% monomial basis
A=vander(x);
% a=T\exp(x')
a=A\(y1(x)'); cond(A)
A2=vander(x2);
y3=A2(:,end-n:end)*(a);

a(end)
%%
clf;
semilogy(abs(a))
% loglog((abs(a2)))
%%
yy1=[y1(x),fliplr(y1(x))]; yy1(n+2)=[]; yy1(end)=[];
a2=fft(yy1)/length(yy1)*2; a2(1)=a2(1)/2; a2(2:2:end)=a2(2:2:end)*-1;
% a2(1:n)-a(1:n)' % 1e-15 for exp so almost the same

clf;
semilogy(abs(a2))
% loglog((abs(a2)));
% loglog((abs(a2)));
hold on;
% loglog([1:n].^(-2));

%%


%
clf;
plot(x2,y1(x2),'-'); hold on;
plot(x2,y2,'-x')
plot(x2,y3,'-o')
%%
clf;
plot(y2-y1(x2))

%%
semilogy(abs(a),'-x')

%%
clf;x2=-1:0.001:1; hold on;
plot(x,cos(37*acos(x)),'x'); 
plot(x2,cos(37*acos(x2)),'-'); 
plot(x2,cos(45*acos(x2)),'-'); 

plot(x,cos(45*acos(x)),'o')
