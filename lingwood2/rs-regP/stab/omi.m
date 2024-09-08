load('startR40.mat'); 
Rar=[40:-1:10];

omiar=[];


for iii=1:length(Rar)

R=Rar(iii);
% clc;
% x=[0.25;-0.05];
% x=[0.2;0.2];
% x=[0.2;0.02];
%  x=[0.15;    0.0];
 x=xa(:,end);
xa=[x];
eps=1e-6;
zh=@(x,y) imag(imagOmega(x,y,bbar,R,data));

%%

for ii=1:10
d2fdx2=(zh(x(1)+eps,x(2))-2*zh(x(1),x(2))+zh(x(1)-eps,x(2)))/eps^2;
d2fdxdy=(zh(x(1)+eps,x(2)+eps)+zh(x(1)-eps,x(2)-eps)-zh(x(1)+eps,x(2)-eps)-zh(x(1)-eps,x(2)+eps))/eps^2/4;
d2fdy2=(zh(x(1),x(2)+eps)-2*zh(x(1),x(2))+zh(x(1),x(2)-eps))/eps^2;

jac=[d2fdx2,d2fdxdy; d2fdxdy, d2fdy2];
g=[(zh(x(1)+eps,x(2))-zh(x(1)-eps,x(2)))/2/eps;(zh(x(1),x(2)+eps)-zh(x(1),x(2)-eps))/2/eps];
x=x-jac\g;
xa=[xa,x];
det(jac);
norm(g);
% imagOmega(x(1),x(2),bbar,R,data)

fprintf("ar: %4.2e\tai: %4.2e\tR: %4.2f\tnorm(g): %4.2e\tomi: %4.2e\n",x(1),x(2),R,norm(g),zh(x(1),x(2)));
if(norm(g)<1e-9), break; end
end
if(norm(g)>1e-9), error("no conv"); end

omiar(iii)=zh(x(1),x(2));
end

%%
plot(Rar,omiar,'-x'); grid on
