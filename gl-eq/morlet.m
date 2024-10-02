addpath     '/people/gesla/Documents/git/rotst2/scripts/source_for_mtb'

x=-3:0.01:3; x=x(1:end-1);
y=sin(5*2*pi*x); y2=exp(-(x*2).^2);

plot(x,y.*y2); 
[f,z]=ft(x,y);
[f2,z2]=ft(x,y2);
[f3,z3]=ft(x,y.*y2);
%%
clf;
plot(f,z,'-x');  xlim([0 9]); grid on; hold on;
plot(f2,z2,'-x');
plot(f3,z3,'-x');
%%
clf;
plot(imag(fft(y.*y2))); xlim([0 50])
