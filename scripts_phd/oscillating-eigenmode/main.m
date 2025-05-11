vr=1;
vi=2;
t=0:0.001:10;
lamr=0;
lami=2;
vp=exp((lamr+lami*i)*t)*(vr+i*vi)+exp((lamr-lami*i)*t)*(vr-i*vi);

%%
plot(t,abs(vp)); hold on; plot(t,2*(cos(lami*t)*vr-vi*sin(lami*t)))

%%
A=[-0.1 1;
    -1 -0.1];
[evc,ev]=eig(A);

dt=0.01;
x=[1;1];
xM=[];

for t=0:dt:100
    x=A*x*dt+x;
    xM=[xM,x];
end

%%
plot(xM(1,:),xM(2,:))