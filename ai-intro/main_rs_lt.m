clc; clear; 
addpath /people/gesla/Documents/git/rotst2/scripts/source_for_mtb


%% 1D

x=0:0.01:2*pi;
nt=3; 
an=rand(nt+1,1)*2-1; an(1)=-sum(an(2:end)*2); bn=rand(nt+1,1)*2-1; bn(1)=0; 
f=an(1)+cos(x'*[1:nt])*an(2:end)*2+sin(x'*[1:nt])*bn(2:end)*2;
plot(x,f)

%% 2D

a=importdata('u20-600-160.dat');
a=loadrs2(a);
x=a.xc; z=a.zc;

%%

ntx=8; ntz=2; lx=10; lz=1;
an=rand(ntz+1,ntx+1)*2-1; %an(1)=-sum(an(2:end)*2); 
bn=rand(ntz+1,ntx+1)*2-1; %bn(1,:)=0; bn(:,1)=0; 

f=0;
% for ikx=0:ntx
%     for ikz=0:ntz
%         % f=f+cos(x/lx*2*pi*ikx).*cos(z/lz*2*pi*ikz)*an(ikz+1,ikx+1);
%         f=f+sin(x/lx*2*pi*ikx).*sin(z/lz*2*pi*ikz)*bn(ikz+1,ikx+1);
%     end
% end

for ikx=-ntx:ntx
    for ikz=-ntz:ntz
        % f=f+cos(x/lx*2*pi*ikx).*cos(z/lz*2*pi*ikz)*an(ikz+1,ikx+1);
        f=f+(an(abs(ikz)+1,abs(ikx)+1)+bn(abs(ikz)+1,abs(ikx)+1)*1i*(1-2*(ikx<0))*(1-2*(ikz<0)))*exp(1i*ikx*x/lx*2*pi).*exp(1i*ikz*z/lz*2*pi);
    end
end


pcolor(x,z,real(f)); shading interp; axis equal; xlim([0 10]); ylim([0 1]); colormap(parula(8)); colorbar();
% mesh(x,z,f); shading interp;  xlim([0 10]); ylim([0 1]); colormap(parula(8));

%harder than I thought

