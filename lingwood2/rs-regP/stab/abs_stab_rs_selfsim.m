addpath ..
clc; clear;
 bbar=0;
% reharr=1050:5:1250;
reharr=200:-2:100;
reharr=1000;
resarr=[reharr(1);0.25;    0.01; 50; 0; 0]; %ig

%%


for ire=1:length(reharr)
    reh=reharr(ire);
a=get_rs(reh,320); data=a;

% clc; 
% x=[0.25;    0.01; 50];
x=resarr(2:4,end);
xa=[x];
eps=1e-6;

zh=@(x,y,R) imag(imagOmega(x,y,bbar,R,data));

%

for inewt=1:15

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
fprintf("reh: %4.2f\t ar: %4.2e\tai: %4.2e\tR: %4.2f\tnorm(g): %4.2e\tomi: %4.2e\n",reh,x(1),x(2),x(3),norm(g),g(end));
if(norm(g)<1e-8) break; end;
end
resarr=[resarr,[reh;x(1);x(2);x(3);norm(g);g(end)]];
if(norm(g)>1e-5) error("no conv"); end

end