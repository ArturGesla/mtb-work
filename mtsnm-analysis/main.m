a=importdata("out-mtsnm");
%%
plot(a(:,1),a(:,12))

%%
ind=10000*3;
x=a(ind:ind+20000,12);
t=a(ind:ind+20000,1);
T=t(end)-t(1);
df=1/T;
z=fft(x-mean(x));
f=0:df:1/mean(diff(t));
semilogy(f(2:end),abs(z(2:end))); hold on;
xlim([0 2])
t(1)

[c,b]=max(abs(z));
f(b)
semilogy(f(b),c,'or')
per=1/f(b)
dt=mean(diff(t))
