% reduced nt version

nt = 7;
nvar = 4;
nr = 70+2;
nz = 100+2;

nv = 4;
ir1 = 36;
jz1 = 68;
ir2 = 36;
jz2 = 15;

ks = -(nt-1):nt-1;

k = 1:2*nt;

time = linspace(0,1,191);
time(end) = [];

%s=read_sunfluidh_probes; plot(s.w(end -100:end,2), s.w(end-100:end,1))


for Rey = 2800:100:4300 % UPO    
   filename = ['XsavRey',num2str(Rey),'BigUPO.mat'];
% for Rey = 2700:100:4400
%     filename = ['XsavRey',num2str(Rey),'Big.mat'];
load(filename)
%load('XsavRey3500Big.mat')
%load('XsavRey2700Big.mat')
%load('XsavfullRe2700.mat')
%load('XsavfullRe2600.mat')

pos1 = nv -1 + (k-1)*nvar + (ir1-1)*nvar*2*nt + (jz1-1)*nvar*nr*2*nt + 1;
pos2 = nv -1 + (k-1)*nvar + (ir2-1)*nvar*2*nt + (jz2-1)*nvar*nr*2*nt + 1;

X1 = X(pos1);
X1 = X1(1:nt) + 1i*X1(nt+1:2*nt);
X1 = [conj(flipud(X1(2:end))) ; X1];

X2 = X(pos2);
X2 = X2(1:nt) + 1i*X2(nt+1:2*nt);
X2 = [conj(flipud(X2(2:end))) ; X2];


xrebuild1 = sum(X1.'.*exp(1j*2*pi*ks.*time'),2);
xrebuild2 = sum(X2.'.*exp(1j*2*pi*ks.*time'),2);

%plot(real(xrebuild2),real(xrebuild1),'k+')
plot3(Rey*ones(size(xrebuild2)),real(xrebuild2),real(xrebuild1))
hold on

end
hold off
%title('Sunfluidh versus Newton Rotor Stator H/R=1.5')
title('PO for Rotor Stator H/R=1.5')

set(gca,'FontSize',16)
%legend('2800','2700','2600')
xlabel('Rey')
ylabel('V_t (r=0.493, z=0.202)') 
zlabel('V_t (r=0.493, z=0.997)') 
