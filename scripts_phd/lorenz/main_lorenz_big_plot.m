r=0:30;
cp=real(sqrt(b*(r-1)));
cm=real(-sqrt(b*(r-1)));

%%
b=8/3;
r2=[1:0.01:24.74];
r1=[0:0.01:1];
r3=[24.74:0.01:30];

close all; hold on;
plot3(r1,0*r1,0*r1,'b-'); hold on;
plot3([r2,r3],0*[r2,r3],0*[r2,r3],'b--','HandleVisibility','off'); hold on;
plot3(r2,real(sqrt(b*([r2]-1))),real(sqrt(b*([r2]-1))),'b-','HandleVisibility','off'); hold on;
plot3(r2,-real(sqrt(b*([r2]-1))),-real(sqrt(b*([r2]-1))),'b-','HandleVisibility','off'); hold on;

plot3(r3,real(sqrt(b*([r3]-1))),real(sqrt(b*([r3]-1))),'b--','HandleVisibility','off'); hold on;
plot3(r3,-real(sqrt(b*([r3]-1))),-real(sqrt(b*([r3]-1))),'b--','HandleVisibility','off'); hold on;

plot3((uMCm(3:neq:end-1,1)*0+1)*rMm(end),uMCm(1:neq:end-1,end),uMCm(2:neq:end-1,end),'r-','HandleVisibility','on');  hold on;
plot3((uMCm(3:neq:end-1,1)*0+1)*rMm(1:11:end),uMCm(1:neq:end-1,1:11:end),uMCm(2:neq:end-1,1:11:end),'r-','HandleVisibility','off');  hold on;

plot3((uMCp(3:neq:end-1,1)*0+1)*rMp(1:20:end),uMCp(1:neq:end-1,1:20:end),uMCp(2:neq:end-1,1:20:end),'r-','HandleVisibility','off');  hold on;

plot3((uMCu(3:neq:end-1,1)*0+1)*rMu(1:10:end),uMCu(1:neq:end-1,1:10:end),uMCu(2:neq:end-1,1:10:end),'g-')

legend("fixed points","upo due to subHopf at r=24.74","upo presented last time by Laurent","Location","best");

grid on; grid minor; xlabel("r"); ylabel("x"); zlabel("y"); title("Biforcation diagram Lorenz b=8/3 sigma=10")