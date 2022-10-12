% a=(vss(180,:))
% close all
% plot(a)
% plot(a>max(a)/2)
b=(abs(a)>(max(a)/2));
c=cumsum(b);
i1=sum(1-cummax(b));
i2=length(r)-sum(1-cummax(b,"reverse"));
di=270;
ic=(i1+i2)/2;
i1=max(round(ic-di/2),1); i2=round(ic+di/2)-min(round(ic-di/2),0);
% hold on;
% plot(r,a)
% plot(r(i1:i2),a(i1:i2)); 