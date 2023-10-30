clear
!rm -rf Usav.mat
disp('Removing USav.mat')

% Finding Periodic Orbits for Langford system with Harmonic Balance Method
%
% We use symbolic computations to derive the Jacobian.
%
% This version split real and imaginary part so that the matrix is real.
% We are thus able to use the fact that Fourier modes are complex
% conjugate and save some memory.
% 
% clear all; close all

%%
%%
close all;

%lbda=1.7; xyz_init = [ 0.2334 -0.9401 1.0858]; Period = 25.05;

lbda = 1.8; xyz_init = [0.76  0.0  1.38197]; Period = 25.13; %8*pi; %25.13;

lbda = 1.97

%% Lorenz 
%     dotX = sigma*(Y - X);
%     dotY = X*(rho - Z) - Y;
%     dotZ = X*Y - beta*Z;
%
% Here we get the UPO with proper initial condition 

% Note : According to Strogatz 
% for sigma = 10  and beta = 8/3,
% Lorenz is Hopf subcritical at rho = 24.74, the UPO branch can be followed
% down to rho =  13.926.

% Classical Lorenz parameters
%sigma = 10;
%beta  = 8/3;

% Some initial guesses
%rho   = 28.; xyz_init =  [15.468238775920206   15.468238775920206   36.546917872586619]; Period=1.5586; 
%rho= 22; xyz_init = [3.029432017323835   3.029432017323835  16.292272216787758]; Period = 1/1.309256812098576;
%rho = 20; xyz_init = [1.633212949507859  1.633212949507859 13.285785046591739]; Period = 1/1.141327248523357;
%rho = 24.; xyz_init =  [5.421095939029748   5.421095939029748  20.323920878698125]; Period = 0.679336764196283;
%rho   = 160.; xyz_init =  [4.1475781218192   5.0524607933205   118.6139607475268]; Period=1.15; 

%rho = 138 change rho for continuation 
% Numerical parameters
ntphys = 100; % number of grid points for initial guess (and plotting)
nt = 12;  % number of modes

itnewtonmax = 8;

%
time = linspace(0,1,ntphys+1);

%tend = 0
%%
options = odeset('RelTol',1e-10);
%Period = 1
%time = linspace(tend,tend+300,ntphys+1);
%time = [tend (tend+300)]
% Integrate in time Lorenz with an initial value on the UPO. 
% We rescaled Lorenz equations so that the orbit has a period = 1.

%lorenz = @(t,x) [Period*sigma*(x(2)-x(1));Period*(rho*x(1)-x(1)*x(3)-x(2));Period*(x(1)*x(2)-beta*x(3))]; % Anonymous Function
%[T,X] = ode45(lorenz, time ,xyz_init); 


langford = @(t,y) [Period*((lbda-3)*y(1)-1/4*y(2)+y(1)*(y(3)+0.2*(1-y(3)^2)));
                   Period*((lbda-3)*y(2)+1/4*y(1)+y(2)*(y(3)+0.2*(1-y(3)^2)));
                   Period*(lbda*y(3)-(y(1)*y(1)+y(2)*y(2)+y(3)*y(3)))        ]; % Anonymous Function
 [T,X] = ode45(langford, time, xyz_init ,options );
%plot(X(:,1),X(:,2))
%plot3(X(:,1),X(:,2),X(:,3)) 

%plot(T,X(:,1))
xyz_init = X(end,:)
tend = T(end)
length(T)
%return


figure(1)
plot3(X(:,1),X(:,2),X(:,3))
hold on;
plot3(X(1,1),X(1,2),X(1,3),'ro')
grid on

figure(2)
plot(T,X)
legend('x','y','z')

figure(3)
plot(T,X(:,3))
hold on
%return
%%

%% Initial guess from ode45.

xguess = X(:,1);
yguess = X(:,2);
zguess = X(:,3);
Fqguess = 1/Period; 

time(end)   = [];   % Remove the last point.
xguess(end) = [];
yguess(end) = [];
zguess(end) = [];

figure(3)
plot(time,xguess,'-')
hold on 

Xf = fft(xguess)/ntphys; % rescale with ntphys 
Yf = fft(yguess)/ntphys;
Zf = fft(zguess)/ntphys;

%return

%% 
% We will pick nt Fourier modes : the mean and the lowest freq in the signal 
ktr = 1:nt;
ks =  0:(nt-1);  

nvar = 3;
ndim = nt*nvar + 1; 

% Declare symbolic unknowns 
syms xpr xpi  ypr ypi zpr zpi [1 nt]
syms Fq  % 1/period
%
% Array of unknowns 
Unknownsv = [xpr ypr zpr xpi ypi zpi]; 
assume(Unknownsv,'real')
% Array of unknowns with Fixed location 
Unknownsf = Fq ;  
assume(Unknownsf,'real')

% Complex variables to simplify equations writing.
xp = xpr + 1i*xpi;
yp = ypr + 1i*ypi;
zp = zpr + 1i*zpi;

Pi = sym(pi);

%%   EQUATIONS 

%% X Equation 
%EqX  = sigma*(yp-xp) - Fq*2*1j*Pi*ks.*xp;  % 
EqX  = (lbda-3)*xp - 0.25*yp + sspconv(xp,zp) - 0.2*sspconv3(xp,zp,zp) + 0.2*xp  - Fq*2*1j*Pi*ks.*xp;
%% Y Equation
%EqY  = rho*xp - sspconv(xp,zp) - yp - Fq*2*1j*Pi*ks.*yp;
EqY  = 0.25*xp + (lbda-3)*yp + sspconv(yp,zp) - 0.2*sspconv3(yp,zp,zp) + 0.2*yp - Fq*2*1j*Pi*ks.*yp;
%% Z Equation
%EqZ  = sspconv(xp,yp) - beta*zp - Fq*2*1j*Pi*ks.*zp;
EqZ  = lbda*zp - (sspconv(xp,xp) + sspconv(yp,yp) + sspconv(zp,zp)) - Fq*2*1j*Pi*ks.*zp;
%% Freq Equation. We impose dotX = 0 at t=0. 
EqFq = sum(ks.*xp) - sum(ks.*conj(xp))  - 0 ;  % This equation is purely imaginary.

%% Nodes location
% X Nodes
Xvarpos = 1;
XNodes   = (0:nt-1) + Xvarpos;
% Y Nodes
Yvarpos = 1+nt;
YNodes   = (0:nt-1) + Yvarpos;
% Z Nodes
Zvarpos = 1+2*nt;
ZNodes   = (0:nt-1) + Zvarpos;
% Freq node
FqNode = ndim;

%Equations= [EqX,EqY,EqZ,EqFq];
Equations= [real(EqX),real(EqY),real(EqZ),imag(EqX),imag(EqY),imag(EqZ),imag(EqFq)];

Unknowns = [Unknownsv,Unknownsf];

Jac = jacobian(Equations,Unknowns);

%return

%% Start Computation
% Initial vector
%Reread initial guess from disk or start from the UPO computed by ode45
if isfile('Usav.mat')
    load('Usav.mat')
    if length(U) ~= ndim
        error('Usav.mat size problem. Try : !rm  Usav.mat') 
    end
    disp(' ')
    disp('*******')
%    disp(['Restarting with rho = ' num2str(rhosave)])
    disp(['Restarting with lbda = ' num2str(lbdasave)])
    disp('*******')
else % use ode45 solution as a guess.
    disp(' ')
    disp('*******')
%    disp(['Starting with guessed solution for rho = ' num2str(rho)])
    disp(['Starting with guessed solution for lbda = ' num2str(lbda)])
    disp('*******')
    U = zeros(ndim,1);
    U(XNodes) = Xf(ktr);  % Get the selected spectral coefficients 
    U(YNodes) = Yf(ktr);  
    U(ZNodes) = Zf(ktr);
    U(FqNode) = Fqguess;
end

% We check by rebuilding in the physical space the solution from the selected
% spectral modes (be carefull between  ' and .')
xrebuild = sum(U(XNodes).'.*exp(1j*2*pi*ks.*time'),2) + sum(U(XNodes)'.*exp( -1j*2*pi*ks.*time'),2)-U(XNodes(1));
plot(time,real(xrebuild))
legend('Solution from ODE',['Rebuild solution with ' num2str(nt) ' modes'])
%return
%%
figure(2)
clf
subplot(311)
xrebuild = sum(U(XNodes).'.*exp(1j*2*pi*ks.*time'),2) + sum(U(XNodes)'.*exp( -1j*2*pi*ks.*time'),2) -U(XNodes(1));
plot(time,real(xrebuild),'b')
xlabel('Initial Guess','FontSize',16)
hold on
subplot(312)
yrebuild = sum(U(YNodes).'.*exp(1j*2*pi*ks.*time'),2) + sum(U(YNodes)'.*exp( -1j*2*pi*ks.*time'),2)-U(YNodes(1));
plot(time,real(yrebuild),'r')
xlabel(['Fq = ' num2str(U(end))],'FontSize',16)
hold on
subplot(313)
zrebuild = sum(U(ZNodes).'.*exp(1j*2*pi*ks.*time'),2) + sum(U(ZNodes)'.*exp( -1j*2*pi*ks.*time'),2)-U(ZNodes(1));
plot(time,real(zrebuild),'g')
hold on

figure(3)
clf
subplot(311)
stem(real(U(XNodes)),'b')
hold on
stem(imag(U(XNodes)),'r+')
xlabel('Initial Guess','FontSize',16)
subplot(312)
stem(real(U(YNodes)),'b')
hold on
stem(imag(U(YNodes)),'r+')
xlabel(['Fq = ' num2str(U(end))],'FontSize',16)
hold on
subplot(313)
stem(real(U(ZNodes)),'b')
hold on
stem(imag(U(ZNodes)),'r+')
hold on
%pause
%return 

%%
for itnewton = 1:itnewtonmax
    
    % Very slow here because no use "matlabfunction" and the matrix is full and no sparsity is assumed !
    tic
    UrUi = ([real(U(1:end-1));imag(U(1:end-1));real(U(end))]);
    Jacd = double(subs(Jac,Unknowns,UrUi.'));
    RHS = -double(subs(Equations,Unknowns,UrUi.')).';
    t=toc;
    disp(['Building matrix ... ',num2str(t),' sec.']);
    %return
    % Solving the system
    tic
    dU = Jacd\RHS;
    t=toc;
    disp(['Solving System ... ',num2str(t),' sec.']);
    
    U = U+[(dU(1:ndim-1)+1i*dU(ndim:2*(ndim-1)));dU(end)];
    
    residualdX = sum(abs(dU)); residualRHS = sum(abs(RHS));
    disp(['                    Sum dX   =   ',num2str(residualdX)])
    disp(['                    Sum RHS  =   ',num2str(residualRHS)])
    
    % Plot solution
    clf
    
    figure(2)
    clf
    subplot(311)
    xrebuild = sum(U(XNodes).'.*exp(1j*2*pi*ks.*time'),2)+ sum(U(XNodes)'.*exp( -1j*2*pi*ks.*time'),2) -U(XNodes(1));
    plot(time,real(xrebuild),'b')
    xlabel(['Iter Newt # = ' num2str(itnewton)],'FontSize',16)
    hold on
    subplot(312)
    yrebuild = sum(U(YNodes).'.*exp(1j*2*pi*ks.*time'),2)+ sum(U(YNodes)'.*exp( -1j*2*pi*ks.*time'),2)-U(YNodes(1));
    plot(time,real(yrebuild),'r')
    xlabel(['Fq = ' num2str(U(end))],'FontSize',16)
    hold on
    subplot(313)
    zrebuild = sum(U(ZNodes).'.*exp(1j*2*pi*ks.*time'),2)+ sum(U(ZNodes)'.*exp( -1j*2*pi*ks.*time'),2)-U(ZNodes(1));
    plot(time,real(zrebuild),'g')
    hold on
    %pause
    %return
    figure(3)
    clf
    subplot(311)
    stem(real(U(XNodes)),'b')
    hold on
    stem(imag(U(XNodes)),'r+')
    xlabel(['Iter Newt # = ' num2str(itnewton)],'FontSize',16)
    subplot(312)
    stem(real(U(YNodes)),'b')
    hold on
    stem(imag(U(YNodes)),'r+')
    xlabel(['Fq = ' num2str(U(end))],'FontSize',16)
    hold on
    subplot(313)
    stem(real(U(ZNodes)),'b')
    hold on
    stem(imag(U(ZNodes)),'r+')
    hold on
    
    
    %pause
    
end % Newton iteration

%% Some extra plotting
disp('Results : ')
disp(['Freq = ' num2str(U(end)) ])

figure(1)
plot3(real(xrebuild),real(yrebuild),real(zrebuild),'r')
set(gca,"FontSize",12,"FontName","Latin Modern Math");
set(gcf,'defaulttextinterpreter','latex')
%title({['Lorenz : $\rho$ = ' num2str(round(rho,4)) ' $\sigma$ = ' num2str(sigma) ...
title({['Langford : $\lambda$ = ' num2str(round(lbda,4)) ...
        ' ntphys = ' num2str(ntphys) ' nt = ' num2str(nt)],  ...
        [' Fq = ' num2str(U(end)) ] ...
        })
%exportgraphics(gcf,'plot.png','Resolution',200)
% axis('equal')    


% Save for next run to use it as guess
%rhosave = rho;
%save('Usav.mat','U','rhosave')
lbdasave = lbda;
save('Usav.mat','U','lbdasave')


%% If needed we can compute the Floquet coefficients.

FloquetCoeff = eig(Jacd(1:end-1,1:end-1));
[~, indices]=sort(abs(imag(FloquetCoeff)));
LeadingFloquet = ((FloquetCoeff(indices(1:3)))).' % Here we can take safely the real part but it is not general
disp(sprintf('%d & %.0f & %.15f & %.15f & %.15f \\\\',nt,lbda,LeadingFloquet))

FloquetCoeff = eig(Jacd(1:end-1,1:end-1));
%[~, indices]=sort(abs(imag(FloquetCoeff)));
%LeadingFloquet = (real(FloquetCoeff(indices(1:3)))).' % Here we can take safely the real part but it is not general
%disp(sprintf('%d & %.3f & %.15f & %.15f & %.15f \\\\',nt,lbda,LeadingFloquet,Periodsave))
radius = 2*abs(U(1+1));
zval   = abs(U((1+2*nt)));
Periodsave = 1/U(end);
disp(sprintf('%d & %.3f & %.15f & %.15f & %.15f \\\\',nt,lbda,radius,zval,2*pi/Periodsave))

save(['Floquetsavlbda' num2str(lbda) 'nt' num2str(nt) '.mat'],'FloquetCoeff','lbdasave','Periodsave')



%% Nested Functions.

function Eqk = sspconv3(u,v,w)
% We assume u and v have the same size and are vectors 
% Here, we return only the 
% central part of the convolution of the same size than u 

% Ugly coding with loops  but OK for now because it is used only in
% preprocessing. Need to clarify the coding 

nt = 2*length(u)-1;
nshift = (nt+1)/2; 
syms Eqk [1 length(u)];
ktab = -(nt-1)/2:(nt-1)/2;

uf = [conj(fliplr(u(2:end))) u ];
vf = [conj(fliplr(v(2:end))) v ];
wf = [conj(fliplr(w(2:end))) w ];

for knum= 0:(nt-1)/2 %ktab
    kind = knum + 1;
    Eqk(kind) = 0;
    for pnum = ktab
        for lnum = ktab
            for qnum = ktab
                if knum==pnum+lnum+qnum
                    Eqk(kind) = Eqk(kind) + uf(pnum+nshift)*vf(lnum+nshift)*wf(qnum+nshift);
                end
            end
        end
    end
end

return
end

function Eqk = sspconv(u,v)
% We assume u and v have the same size and are vectors 
% Here, we return only the 
% central part of the convolution of the same size than u 

% Ugly coding with loops  but OK for now because it is used only in
% preprocessing. Need to clarify the coding 

nt = 2*length(u)-1;
nshift = (nt+1)/2; 
syms Eqk [1 length(u)];
ktab = -(nt-1)/2:(nt-1)/2;

uf = [conj(fliplr(u(2:end))) u ];
vf = [conj(fliplr(v(2:end))) v ];

for knum= 0:(nt-1)/2 %ktab
    kind = knum + 1;
    Eqk(kind) = 0;
    for pnum = ktab
        for lnum = ktab
            if knum==pnum+lnum
                Eqk(kind) = Eqk(kind) + uf(pnum+nshift)*vf(lnum+nshift);
            end
        end
    end
end

return
end

