clear all;
sigma=10; r=21; b=8/3;
%tr ch 13.928 to 24.06
lorenz = @(t,y) [sigma*(y(2)-y(1));
    y(1)*(r-y(3))-y(2);
    y(1)*y(2)-b*y(3)];

nt=10000; dt=0.01; t0=0;
% x0=[-2,0,0]';[t,y]=rk4_2(@(t,y)lorenz(t,y),nt,dt,t0,x0);

c1=[sqrt(b*(r-1));sqrt(b*(r-1));r-1];
c2=[-sqrt(b*(r-1));-sqrt(b*(r-1));r-1];

% en=min(vecnorm(y-c2,1),vecnorm(y-c1,1));

np=2000;
enArr=zeros(np,nt+1);

rng(1);
x0a=[];
for i=1:np
%     x0=c1+(rand(3,1)-0.5)*20;
%     x0=c1+([0;0;30+rand(1,1)-0.5]);
%     x0=c1+([0;0;rand(1,1)-0.5]*12+[0;0;14]);
    x0=[c1(1:2);32+2/np*i];
    x0a=[x0a;x0'];
    [t,y]=rk4_2(@(t,y)lorenz(t,y),nt,dt,t0,x0);
    en=min(vecnorm(y-c2,1),vecnorm(y-c1,1));
%     enArr=[enArr;en];
    enArr(i,:)=en;
    i
end
%%
% plot(t,y')
% semilogy(t,en)
clf;
ev=100;
plot(t(1:ev:end),enArr(1:ev:end,1:ev:end)'); grid ;
% semilogx(t,enArr');
% semilogy(t,enArr')
xlabel("time"); ylabel("distance from a laminar point"); title("Lorenz r=21 | Lifetimes with different z pert")
exportgraphics(gcf,"ltlor1.png")

%%
b=enArr>10;
lt=sum(b')*dt;

lts=sort(lt);

bb=fliplr(lts);
bb2=1:length(lts);

semilogy(bb,bb2/length(bb2),'-x')
xlabel("time"); ylabel("P(lifetime > time)"); title("Lorenz r=21 | Lifetimes with different z pert")
exportgraphics(gcf,"ltlor2.png")

%%

plot3(x0a(:,1),x0a(:,2),x0a(:,3),'x'); hold on;
plot3(c1(1),c1(2),c1(3),'r+');
plot3(c2(1),c2(2),c2(3),'r+');

%%
dlmwrite('lorenzdata.dat',[x0a,lt'],'delimiter',' ');
%%
close all;
[a,b]=sort(x0a(:,3));
plot(x0a(b,3),lt(b),'.-');
xlabel("z0"); ylabel("lifetime"); title("Lorenz r=21 | Lifetimes with different z pert")
exportgraphics(gcf,"ltlor3.png")
