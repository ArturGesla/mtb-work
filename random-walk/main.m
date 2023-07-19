clear; clc; close all;
% hold on;
% for N=2:4:200
for N=10:10:40
    y=[];
    xA=[];
    for x=-N:2:N
        xA(end+1)=x;
        y(end+1)=prob(N,x);
    end
    plot(xA,y); hold on; grid on;
end