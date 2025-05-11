clc; clear; 
load('neutcurve200.mat');


% m=0:32;
% bbarar=[0:0.01:0.7];
bbarar=[0:-0.01:-0.7];
recar=[];

%%
for ib=1:length(bbarar)

x=xa(:,end);
xa=[x];
bbar=bbarar(ib);


eps=1e-6;
zh=@(x,y,R) imag(imagOmega(x,y,bbar,R,data));

for ii=1:10
d2fdx2=(zh(x(1)+eps,x(2),x(3))-2*zh(x(1),x(2),x(3))+zh(x(1)-eps,x(2),x(3)))/eps^2;
d2fdxdy=(zh(x(1)+eps,x(2)+eps,x(3))+zh(x(1)-eps,x(2)-eps,x(3))-zh(x(1)+eps,x(2)-eps,x(3))-zh(x(1)-eps,x(2)+eps,x(3)))/eps^2/4;
d2fdy2=(zh(x(1),x(2)+eps,x(3))-2*zh(x(1),x(2),x(3))+zh(x(1),x(2)-eps,x(3)))/eps^2;

d2fdxdb=(zh(x(1)+eps,x(2),x(3)+eps)+zh(x(1)-eps,x(2),x(3)-eps)-zh(x(1)+eps,x(2),x(3)-eps)-zh(x(1)-eps,x(2),x(3)+eps))/eps^2/4;
d2fdydb=(zh(x(1),x(2)+eps,x(3)+eps)+zh(x(1),x(2)-eps,x(3)-eps)-zh(x(1),x(2)+eps,x(3)-eps)-zh(x(1),x(2)-eps,x(3)+eps))/eps^2/4;

dfdx=(zh(x(1)+eps,x(2),x(3))-zh(x(1)-eps,x(2),x(3)))/eps/2;
dfdy=(zh(x(1),x(2)+eps,x(3))-zh(x(1),x(2)-eps,x(3)))/eps/2;
dfdb=(zh(x(1),x(2),x(3)+eps)-zh(x(1),x(2),x(3)-eps))/eps/2;

jac=[d2fdx2,d2fdxdy,d2fdxdb;
    d2fdxdy, d2fdy2,d2fdydb;
    dfdx,dfdy,dfdb];

% g=[(zh(x(1)+eps,x(2))-zh(x(1)-eps,x(2)))/2/eps;(zh(x(1),x(2)+eps)-zh(x(1),x(2)-eps))/2/eps];

g=[dfdx;dfdy;zh(x(1),x(2),x(3))];

x=x-jac\g;
xa=[xa,x];
det(jac);
norm(g);
fprintf("bbar: %4.2e\tar: %4.2e\tai: %4.2e\tR: %4.2f\tnorm(g): %4.2e\tomi: %4.2e\n",bbar,x(1),x(2),x(3),norm(g),g(end));

if(norm(g)<1e-9), break; end
end
if(norm(g)>1e-9), error("no conv"); end

recar(ib)=x(end);
fprintf("\n");
end

%%
% save("betaplusreh200.mat",'bbarar','recar')
save("betaminusreh200.mat",'bbarar','recar')
%%

plot(recar,bbarar(1:length(recar)));
%%

 addpath C:\Users\Artur\Documents\GitHub\rotst2\scripts\source_for_mtb;
 
 clf; hold on;
a=load("betaminusreh200.mat");   plot(a.recar,a.bbarar(1:length(a.recar)),'-k');
a=load("betaplusreh200.mat");   plot(a.recar,a.bbarar(1:length(a.recar)),'-k');

grid on; title("Absolute inst. threshold Reh=200"); xlim([50,150]); ylim([-0.1 0.25]);

for m=-40:2:40
    rr=20:200;
    plot(rr,m./rr,'-',"Color",[0.7 0.7 0.7]);    
end

for m=-12:4:32 
    text(140,m/140,"m="+num2str(m))
end
fnts=12; xlabel("$Re_c$"); ylabel("$\beta=m/Re$"); jfm_plt_aid_comm;