clc; clear;
t=linspace(-1,1,100);

% plot(x,y);

% y=@(x) sin(x*3/2*pi)+0.2*cos(x*3/2*pi);
y=@(x) cos(x*3/2*pi);
% y=@(x) x.^2;
% y=@(x) x.^4;
% y=@(x) x+ 2*x.^2-1+1;
yt=y(t);
%
cnt=8;
tch=cos(linspace(0,pi,cnt));
yc=y(tch);

yc2=[yc,fliplr(yc(2:end-1))];
z=fft(yc2)/length(yc2);

ych=0;
for i=1:cnt
    ych=ych+real(z(i))*(1+(i~=1)*1)*cos((i-1)*acos(t));
end


%
clf; close all;
plot(t,yt); hold on;
plot(tch,yc,'x');
plot(t,ych,'--');
save("my.mat")

%% laurent mtd of interpolation

clear;
cnt=8;
tch=cos(linspace(0,pi,cnt));

y=@(x) sin(x*3/2*pi)+0.2*cos(x*3/2*pi);
y=@(x) cos(x*3/2*pi);

% y=@(x) x+ 2*x.^2-1+1;

A=cos(acos(tch)'*[0:cnt-1]);
b=y(tch');
z=A\b;

ych=0;
t=linspace(-1,1,100);
for i=1:cnt
    ych=ych+real(z(i))*cos((i-1)*acos(t));
end


%
clf; close all;
plot(t,y(t)); hold on;
plot(tch,y(tch),'x');
plot(t,ych,'--');
save("lau.mat")

%%
my=load("my.mat");

lau=load("lau.mat");
diff([my.z(2:8)'*2, lau.z(2:end)]')
