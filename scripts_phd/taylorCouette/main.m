% i=430;
cd     'C:\Users\Artur\Documents\taylorCouette'
%%
i=2610

clc; clf;
% for i=1500
a=importdata("./figs/p."+num2str(i,'%05d')+".bmp");
% a=importdata("./figs/p1.01200.bmp");
addpath(    'C:\Users\Artur\Documents\GitHub\rotst2\scripts\source_for_mtb')
ix=400; iy=50;
% pcolor((a(1:200,:,1)));shading interp; hold on; plot(ix,iy,'+')
% plot(a(200,200:800,1))
% imshow(a(1:200,:,1))
% imshow(a(1100:1500,360:800,1))
% imshow(a(800:1700,200:950,1))
da=(a(:,:,2)); %whole
% da=(a(600:1800,100:900,3)); %bulk 
% da=(a(1200:1300,500:600,3)); %detail
% da=(a(1200:1200+64,500:500+64,3)); %detail2
% da=(a(1060:1260,640:740,1)); %detail3
% da=(a(1:200,200:950,2)); %top
% da=((a(800:1700,200:950,1))+(a(800:1700,200:950,2))+(a(800:1700,200:950,3)))/3;
% da=flipud(da);
h=pcolor(da); h.EdgeAlpha=0;
colormap(gray(256)); caxis([0 255]); %colorbar();
% shading interp;
hold on; plot(ix,iy,'+')
% imshow(da)
i=i+1;
% ylim([1 300])
title(num2str(i))
% exportgraphics(gcf,"./fig2/p2."+num2str(i,'%05d')+".png");
i
% end

%%
arr=[];
tarr=[];
di=100;
dt=1/30;
% for i=500+di:1100+di
    for i=1:2610
    
    a=importdata("./figs/p."+num2str(i,'%05d')+".bmp");
arr=[arr; a(iy,ix,1)];
tarr=[tarr; i*dt];
i
end

%%
clf;
plot(tarr,arr)
%%
loglog(abs(fft(arr)))
%%
[f,z]=ft(tarr,arr);
%%
semilogy(1./f,z)

%% conv

t=-5:0.1:5;

x=exp(-t.^2);
x2=exp(-(t-2).^2);

%% fft by eye
tt=[8 23 37 53 70];
Tt=[1 1.4 2 3.3 7];

%%

% T=2.3:0.1:10;
T=Tt;
R1=60/2*1e-3;
R2=75/2*1e-3;

% R1=14*1e-3;
% R2=33*1e-3;
% T=2*pi/0.086;
% T=10;

mu=0;
nu=1e-6;
eta=R1/R2;
om1=2*pi./T;
Ta=4*om1.^2*R1^4/nu^2*(1-mu)*(1-mu/eta^2)/(1-eta^2)^2;
format shortG
tt
Re1=om1*R1*(R2-R1)/nu


%%
close all;
x=-5:0.1:5;
da=exp(-x.^2-x'.^2);
da1=exp(-(x-2).^2-x'.^2);
%%
 figure(); imshow(da);
 figure(); imshow(da1);
 DA=da;
 DA1=da1;
 %%
 c=conv2(da-mean(da,'all'),da1-mean(da1,'all')); [e,f]=max(c,[],'all','linear'); ii=mod(f,length(c)); jj=(f-ii)/length(c)+1;
 djj=jj-length(da); dii=ii-length(da);
 %%
 iia=-length(da)+1:1:length(da)-1;
 mesh(iia,iia,c)

 %% piv
 
 winsize=31; halfsize=(winsize+1)/2; hm=(winsize-1)/2;
 
 arr=[];
 
 sz1=length(halfsize:size(DA,1)-hm);
 sz2=length(halfsize:size(DA,2)-hm);
 
 for ix=halfsize:size(DA,1)-hm
     for iy=halfsize:size(DA,2)-hm
     da=DA(ix-hm:ix+hm,iy-hm:iy+hm);
     da1=DA(ix-hm:ix+hm,iy-hm:iy+hm);
     
      c=conv2(da-mean(da,'all'),da1-mean(da1,'all')); [e,f]=max(c,[],'all','linear'); ii=mod(f,length(c)); jj=(f-ii)/length(c)+1;
 djj=jj-length(da); dii=ii-length(da);
 
 arr=[arr; ix,iy,dii,djj;];
 
 
     end
     ix
 end
 
 X=reshape(arr(:,1),[sz2,sz1]);
 Y=reshape(arr(:,2),[sz2,sz1]);
 U=reshape(arr(:,3),[sz2,sz1]);
 V=reshape(arr(:,4),[sz2,sz1]);
 
 %%
 close all;
 quiver(X,Y,U,V);
%  pcolor(X,Y,sqrt(U.^2+V.^2)); shading interp;
     
 
