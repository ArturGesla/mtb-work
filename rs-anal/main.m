clc; clear;

zj=[1	2.4048	3.8317	5.1356	6.3802	7.5883	8.7715;
2	5.5201	7.0156	8.4172	9.7610	11.0647	12.3386;
3	8.6537	10.1735	11.6198	13.0152	14.3725	15.7002;
4	11.7915	13.3237	14.7960	16.2235	17.6160	18.9801;
5	14.9309	16.4706	17.9598	19.4094	20.8269	22.2178];
zj=zeros(100,7);
zj(:,3)=besselzeroj(1,100);

besselj(1,zj(5,3));

L=1;
R=10;
z=0:0.01/2:L;
x=0:0.1:R;
phi=zeros(length(z),length(x));

% A1=@(n)2/(R^2*besselj(2,zj(n,3))*sinh(L/R*zj(n,3)))*(-1/4*R^4*hypergeom([],3,-R^2/4) / gamma(3) );
A1=@(n)2/(R^2*besselj(2,zj(n,3))^2*sinh(L/R*zj(n,3)))*(R^3*besselj(2,zj(n,3))/zj(n,3));
B1=@(n)R*(1-(-1)^n)*2/n/pi/besseli(1,n*pi/L*R,1);
for ix=1:length(x)
    for iz=1:length(z)
        for in=1:15
            phi(iz,ix)=phi(iz,ix)+sinh(z(iz)/R*zj(in,3))*besselj(1,zj(in,3)*x(ix)/R)*A1(in);
            kn=in*pi/L;
            phi(iz,ix)=phi(iz,ix)+sin(kn*z(iz))*besseli(1,kn*x(ix),1)*B1(in)*exp(kn*(x(ix)-R));
        end
    end
end
    %
    close all;mesh(phi)
