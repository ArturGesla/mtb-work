n=15;
x=-cos((0:n)/n*pi);
% y1=x;
y1=@(x) exp(x);
a=[];
T=cos([0:n]'*acos(x));
% a=T\exp(x')
a=T'\(y1(x)')

%%

x2=-1:0.01:1;
y2=a'*cos([0:n]'*acos(x2))
%
clf;
plot(x2,y1(x2),'-'); hold on;
plot(x2,y2)
%%
clf;
plot(y2-y1(x2))

%%
semilogy(abs(a),'-x')