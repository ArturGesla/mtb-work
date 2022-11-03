mu=0.005;
omega=1;
gamma=-10;
tau=0.20; %2*pi/((omega-gamma*mu));
%
J=[mu, -omega; omega,mu];
k=0.045;
% tau=1; 
beta=0;
B=[1+k*tau*cos(beta),-k*tau*sin(beta);k*tau*sin(beta),1+k*tau*cos(beta)];
eig(J,B)

%% analitical
beta=[0:0.1:2*pi]; beta=0;
k=-8:0.1:800;
a=1+k.*tau.*cos(beta);
b=k.*tau.*sin(beta);
delta=(-2*a*mu-2*b*omega).^2-4*(a.^2+b.^2).*(omega^2+mu^2);
sol1=-((-2*a*mu-2*b*omega)-sqrt(delta))/2./(a.^2+b.^2)
sol2=-((-2*a*mu-2*b*omega)+sqrt(delta))/2./(a.^2+b.^2)

plot(k,real(sol1));