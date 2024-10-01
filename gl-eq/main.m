l=1; n=50; dx=l/n; N=n+2; xc=-dx/2:dx:l+dx/2; 
U=1; D=0.1; Pe=U*dx/D;
%%
D1=(diag(-1/2/dx*ones(N-1,1),-1)+diag(1/2/dx*ones(N-1,1),1)); 
D2=D*(diag(1/dx/dx*ones(N-1,1),-1)-2*diag(1/dx/dx*ones(N,1),0)+diag(1/dx/dx*ones(N-1,1),1)); 

A=-U*D1+D2;
% A([1,N],:)=A([1,N],:)*0; A(1,1:2)=[1,1]; A(end,end-1:end)=[1,1]; %dir
% A([1,N],:)=A([1,N],:)*0; A(1,1:2)=[-1,1]; A(end,end-1:end)=[-1,1]; %Neu
A([1,N],:)=A([1,N],:)*0; A(1,1)=1; A(1,end-1)=-1;  A(end,2)=1; A(end,end)=-1;  %per
b=zeros(N,1); b(end)=2;
%% ss
u=A\b; xc=-dx/2:dx:l+dx/2; 
plot(xc,u,'-x'); grid on;

%% ti
u0=rand(N,1);
u=u0;
nt=1000;
dt=0.001;
CFL=U*dt/dx; Dn=D*dt/dx/dx;
uarr=[];
%
for it=1:nt
    f=A*u;
    u(2:end-1)=u(2:end-1)+dt*f(2:end-1); u(1)=u(end-1); u(end)=u(2); %EE
    uarr=[uarr,u];
    norm(u)
end
    
 %%
 plot(xc,u,'-x'); grid on;
 
 %%
mesh(uarr)