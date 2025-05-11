
bbara=0.1:0.02/2:0.18;
Ra=510:10/2:560; Ra=flip(Ra);

[X,Y]=meshgrid(bbara,Ra);

Y(:,2:2:end)=flipud(Y(:,2:2:end));

x=reshape(X,[length(bbara)*length(Ra),1]);
y=reshape(Y,[length(bbara)*length(Ra),1]);

plot(x,y,'-x')

z=x.^2+y.^2*1e-6;
plot3(x,y,z);

X2=reshape(x,[length(Ra),length(bbara)]);
Y2=reshape(y,[length(Ra),length(bbara)]);
Z2=reshape(z,[length(Ra),length(bbara)]);

Y2(:,2:2:end)=flipud(Y2(:,2:2:end));
Z2(:,2:2:end)=flipud(Z2(:,2:2:end));

%%
contour(X2,Y2,Z2)

%%
z=z*0; ar0=z; ai0=z; omr0=z; ng=z;
save('xybbarRa.mat','x','y','z','ar0','ai0','omr0','ng');

%%

is=1; %global iterationindex on the curve
%%
is=is+1;
%%
% for is=5:8
    for is=81:99
%
a=load('xybbarRa.mat');
bbar=a.x(is); R=a.y(is);
fprintf("is: %d \t bbar: %f\t R: %f\n",is,bbar,R);

%
% clc;
% x=[0.25;-0.05];
% x=[0.2;0.2];
x=[0.2;-0.1];
if(is~=1), x=[a.ar0(is-1);a.ai0(is-1);]; end
xa=[x];
eps=1e-6;
zh=@(x,y) imagOmega(x,y,bbar,R);

%
for i=1:10

    zp=zh(x(1),x(2));
d2fdx2=(zh(x(1)+eps,x(2))-2*zp+zh(x(1)-eps,x(2)))/eps^2;
d2fdxdy=(zh(x(1)+eps,x(2)+eps)+zh(x(1)-eps,x(2)-eps)-zh(x(1)+eps,x(2)-eps)-zh(x(1)-eps,x(2)+eps))/eps^2/4;
d2fdy2=(zh(x(1),x(2)+eps)-2*zp+zh(x(1),x(2)-eps))/eps^2;

jac=[d2fdx2,d2fdxdy; d2fdxdy, d2fdy2];
g=[(zh(x(1)+eps,x(2))-zh(x(1)-eps,x(2)))/2/eps;(zh(x(1),x(2)+eps)-zh(x(1),x(2)-eps))/2/eps];
x=x-jac\g;
xa=[xa,x];
% det(jac)
% norm(g)

fprintf("i: %d ng: %1.3e det jac: %1.3e omega0i: %1.3e\n",i,norm(g),det(jac),zp);

if(norm(g)<1e-9), break; end

end

if(norm(g)>1e-5), error("no conv"); end
%
z=a.z; z(is)=zh(x(1),x(2)); 
ar0=a.ar0; ar0(is)=xa(1,end);
ai0=a.ai0; ai0(is)=xa(2,end);
omr0=a.omr0; 
ng=a.ng; ng(is)=norm(g);

x=a.x; y=a.y; 
save('xybbarRa.mat','x','y','z','ar0','ai0','omr0','ng');


end
%%

a=load('xybbarRa.mat');
% is=11;
plot3(a.y(1:is),a.x(1:is),a.z(1:is),'-x');
%%
X2=reshape(a.x,[length(Ra),length(bbara)]);
Y2=reshape(a.y,[length(Ra),length(bbara)]);
Z2=reshape(a.z,[length(Ra),length(bbara)]);

Y2(:,2:2:end)=flipud(Y2(:,2:2:end));
Z2(:,2:2:end)=flipud(Z2(:,2:2:end));

%%
 
contour(Y2,X2,Z2,20,'--'); hold on; colormap(hsv(8)); colorbar();
contour(Y2,X2,Z2,[0 0],'-k'); legend("$\omega_{0,i}$","$\omega_{0,i}=0$")

xlim([500 560]); title("Absolute stability map | $Re_c\approx510$")

xlabel("Re"); ylabel("$\bar{\beta}$")
fnts=10; jfm_plt_aid_comm;

%%
exportgraphics(gcf,"absinst.eps")