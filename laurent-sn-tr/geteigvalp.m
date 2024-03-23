Ray = [2.4e6:1.e5:2.6e6 2.8e6:2e5:5.8e6 6.0e6:0.5e6:15.e6]; % B1
%Ray = 2.4e6; % B1
nx = 128;
ny = 128;
eigB1 = [];
Ray1 = [];
eigB2 = [];
Ray2 = [];
eigB3 = [];
Ray3 = [];
FreqB1 = [];
FreqB2 = [];
NormTB1 = [];
NormTB2 = [];

%Read a first time to get geometry and problem 
load(['EigRay',num2str(Ray(1)),'_',num2str(nx),'x',num2str(ny),'.mat'])
matvol = diff(geo.xu).*diff(geo.yv).';

figure(1)
for Ra = Ray
    name = ['EigRay',num2str(Ra),'_',num2str(nx),'x',num2str(ny),'.mat']
    load(name)
    plot(lbda,'+')
    hold on
    eigB1 = [eigB1,  lbda( (real(lbda)>-0.01 & real(lbda)<0.1874 & imag(lbda)>1.16 & imag(lbda)<1.3311)  ) ];
    Ray1 = [Ray1, Ra];
    eigB2 = [eigB2,  lbda( (real(lbda)>-0.01 & real(lbda)<0.1874 & imag(lbda)>1.3380 & imag(lbda)<1.5261)  ) ];
    Ray2 = [Ray2, Ra];
    B3ext = [lbda(  (real(lbda)>-0.024 & real(lbda)<=0.081 & imag(lbda)>1.642 & imag(lbda)<1.72) ...
                   |  (real(lbda)>0.080 & real(lbda)<=0.21 & imag(lbda)>1.51 & imag(lbda)<=1.642)  ) ];    
    eigB3 = [eigB3, B3ext];
    Ray3 = [Ray3, Ra];

    load(strcat('XsavRay',num2str(Ra),'_',num2str(geo.nx-2),'x',num2str(geo.ny-2),'.mat'));
    Xbase = X;
    ndim = geo.nx*geo.ny*prob.nvar*2;
    f.Wb =  reshape(Xbase(4:2*prob.nvar:end),geo.nx,geo.ny);
% 
%     if prob.nt>1
%         ndim = nx*ny*nvar*2*nt+1;   % The '2' comes from nt real + nt imag
%         nbRHS = 1;
%     else
%         ndim = nx*ny*nvar*2;  % It includes an imag part for consistancy
%         nbRHS = 1;
%     end



    nameB1 = ['B1/XsavRay',num2str(Ra),'_',num2str(nx),'x',num2str(ny),'BigUPO.mat']
    load(nameB1)
    ndim = length(X);  
    nt = (ndim - 1)/ (geo.nx*geo.ny*prob.nvar*2);
    FreqB1 = [FreqB1, X(end)];
    Xinter = reshape(X(1:ndim-mod(1,nt)),prob.nvar,[]);
    Xm0 = reshape(Xinter(:,1:2*nt:end),prob.nvar*geo.nx*geo.ny,1);
    f.Wb1 = reshape(Xm0(4:prob.nvar:end),geo.nx,geo.ny);
    NormTB1 =[NormTB1, sum(sum(abs(f.Wb1(2:end-1,2:end-1) - f.Wb(2:end-1,2:end-1)).*matvol ))];



    nameB2 = ['B2/XsavRay',num2str(Ra),'_',num2str(nx),'x',num2str(ny),'BigUPO.mat']
    load(nameB2)
    ndim = length(X);  
    nt = (ndim - 1)/ (geo.nx*geo.ny*prob.nvar*2);
    FreqB2 = [FreqB2, X(end)];
    Xinter = reshape(X(1:ndim-mod(1,nt)),prob.nvar,[]);
    Xm0 = reshape(Xinter(:,1:2*nt:end),prob.nvar*geo.nx*geo.ny,1);
    f.Wb2 = reshape(Xm0(4:prob.nvar:end),geo.nx,geo.ny);
    NormTB2 =[NormTB2, sum(sum( abs(f.Wb2(2:end-1,2:end-1) - f.Wb(2:end-1,2:end-1)).*matvol ))];


    

end


plot(eigB1,'go')
plot(eigB2,'ro')
plot(eigB3,'bo')
grid on
%%
%  figure(2)
%  plot(Ray1,imag(eigB1)/2/pi,'r+')
%  hold on
%  plot(Ray,FreqB1,'r*')
%  plot(Ray2,imag(eigB2)/2/pi,'b+')
%  plot(Ray,FreqB2,'b*')
%  plot(Ray3,imag(eigB3)/2/pi,'g+')
%  grid on
% %%
%  figure(3)
%  plot(Ray1,-imag(eigB1)/2/pi+FreqB1,'r+')
%  hold on
%  plot(Ray2,-imag(eigB2)/2/pi+FreqB2,'b+')
%  legend('B1','B2')
%  xlabel('Ra')
%  ylabel('Freq HBM - Freq LSA ')
%  title(['Grid cosine ',num2str(nx),'x',num2str(ny)] )

%%
%%


 figure(4)
 plot(Ray1,NormTB1,'r+')
 hold on
 plot(Ray2,NormTB2,'b+')
 legend('B1','B2')
 xlabel('Ra')
 ylabel('Norm (T-Tbase)  ')
 title(['Grid cosine ',num2str(nx),'x',num2str(ny)] )



