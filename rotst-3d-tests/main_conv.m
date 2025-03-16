x=[1,2,3,0,0];
y=[4,5,6,0,0];
conv(x,y)
%%
% fftshift
(ifft(fft(x).*fft(y)))
%works but gl to gen it