clc; clear;
cd     '/people/gesla/Documents/git/mtb-work/langford-system/laurent-sn'

load("Floquetsavlbda1.97nt2.mat"); fl2=FloquetCoeff;
load("Floquetsavlbda1.97nt3.mat"); fl3=FloquetCoeff;
load("Floquetsavlbda1.97nt4.mat"); fl4=FloquetCoeff;


load("eve.mat"); fl1=evE;


cd     '/people/gesla/Documents/git/mtb-work/langford-system/fd_newt'

load("exp-fd-nt-11.mat"); fl11=exp;
load("exp-fd-nt-21.mat"); fl21=exp;
load("exp-fd-nt-31.mat"); fl31=exp;
load("exp-fd-nt-51.mat"); fl51=exp;
load("exp-fd-nt-101.mat"); fl101=exp;
load("exp-fd-nt-201.mat"); fl201=exp;
load("exp-fd-nt-401.mat"); fl401=exp;
%%
close all;

plot(real(fl1),imag(fl1),'o'); grid on; hold on;
% plot(real(fl2),imag(fl2),'x'); grid on; hold on;
% plot(real(fl3),imag(fl3),'s'); grid on; hold on;
% plot(real(fl4),imag(fl4),'>'); grid on; hold on;


% plot(real(fl11),imag(fl11),'+'); grid on; hold on;
% plot(real(fl21),imag(fl21),'+'); grid on; hold on;
% plot(real(fl31),imag(fl31),'+'); grid on; hold on;

plot(real(fl51),imag(fl51),'+'); grid on; hold on;
plot(real(fl101),imag(fl101),'+'); grid on; hold on;
plot(real(fl201),imag(fl201),'+'); grid on; hold on;
plot(real(fl401),imag(fl401),'+'); grid on; hold on;


