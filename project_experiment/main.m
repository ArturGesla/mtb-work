clc; clear;

R=86e-3/2; %containder radius
nuarr=[];
waterMlArr=[];

% glycerineM3=pi*R^2*10e-3;
% glycerineMl=glycerineM3*1000*1000;
glycerineMl=200;


for i=0:1:67
waterMl=i;
volumeRatio=glycerineMl/(waterMl+glycerineMl);

[rho,eta]=density_viscosity_glycerine_mix(volumeRatio,20);
nu=eta/rho;

nuarr=[nuarr;nu];
waterMlArr=[waterMlArr;waterMl];
end

%%

% semilogy(waterMlArr,nuarr)
semilogy(glycerineMl./(waterMlArr+glycerineMl),nuarr,'-x');
grid on;
xlabel("glycerine volume ratio"); ylabel("kinematic viscosity");
exportgraphics(gcf,"vesc.png")
%%
Omega=115*2*pi;
Re=Omega*R^2./nuarr;

plot(waterMlArr,Re);

%
H=10:10:30; H=H*1e-3;
G=H./R;
[RE,GG]=meshgrid(Re,G);
plot(RE,GG,'x');

%
[wml,hh]=meshgrid(waterMlArr,H/1e-3);
plot(wml,hh,'kx'); grid on;
xlabel("water [ml]"); ylabel("cavity height [mm]");
ll=size(wml,1)*size(wml,2);
text(reshape(wml,[ll,1]),reshape(hh,[ll,1]),num2str([1:ll]'))
%%
[waterMlArr,glycerineMl./(waterMlArr+glycerineMl),nuarr,Re]