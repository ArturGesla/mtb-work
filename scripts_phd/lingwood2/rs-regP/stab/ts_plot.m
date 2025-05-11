%%

set(gcf,"Position",[          24         267        1306         399]);
a1=load("resarrto148.mat"); a=a1.resarr(:,2:end);
a1=load("resarrto2000.mat"); a=[fliplr(a),a1.resarr(:,2:end)];

subplot(1,4,1:2); plot(a(1,:),a(4,:),'-'); title("Red vs Reh at abs inst");
xlabel("Reh"); ylabel("Red"); 
hold on;
plot(a(1,:),10*sqrt(a(1,:)),'-');
plot(a(1,:),a(1,:)./a(1,:)*24/sqrt(0.313),':k');
plot(1000,48.5,'r+');
plot(1000,45.53,'bx');


legend("Red crit","Red edge of cavity","Bodewadt limit (LW97)","Serre et al. 2004","mesh double","Location","best"); grid on;

axes('Position',[.17 .6 .05 .25]); 
plot(a(1,:),a(4,:),'-'); hold on; xlim([130 200]); ylim([0 150]);
plot(a(1,:),10*sqrt(a(1,:)),'-'); box on; grid on; 

axes('Position',[.25 .27 .05 .25]); 
plot(a(1,:),a(4,:),'-'); hold on; xlim([900 1100]); ylim([40 50]);
plot(a(1,:),10*sqrt(a(1,:)),'-'); box on; grid on; 
plot(a(1,:),a(1,:)./a(1,:)*24/sqrt(0.313),':k');
plot(1000,48.5,'r+');
plot(1000,45.53,'bx');


subplot(1,4,3); plot(a(1,:),a(2,:)); title("ar vs Reh at abs inst");
xlabel("Reh"); ylabel("ar"); grid on;

subplot(1,4,4); plot(a(1,:),a(3,:)); title("ai vs Reh at abs inst");
xlabel("Reh"); ylabel("ai"); grid on;


%% 3 curve

clf;
cd     'C:\Users\Artur\Documents\GitHub\mtb-work\lingwood2\rs-regP\stab'

a1=load("resarrto148.mat"); a=a1.resarr(:,2:end);
a1=load("resarrto2000.mat"); a=[fliplr(a),a1.resarr(:,2:end)];

 hold on;
r=0:0.01:10; rlw97=24/sqrt(0.313); reh2=(rlw97./r).^2; plot(r,reh2);
red=a(4,:); reh=a(1,:); plot(red./sqrt(reh),reh,'-');

%ti
cd       'C:\Users\Artur\Documents\rotst\self-sim'
R=[]; ress=[];
a=load("dataR10.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR8.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR7.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR6.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR5.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR4.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR3.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR2.5.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR2.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR1.7.mat"); R=[R;a.R]; ress=[ress;a.re0];




plot(R,ress,'-x');
 set(gca,"YTick",[(rlw97/10).^2,100:100:600,900:300:3000]); grid on; hold on;
% set(gca,"XTick",sort([get(gca,"XTick"),1.7,0.78]))
set(gca,"XTick",sort([3:1:10,2.5,1.7,0.78]))

ylabel("Reh"); xlabel("r"); xlim([0 10]); ylim([10 3e3]);
set(gca,"yscale","log");
legend("assuming Bodewadt profile","taking the real ss profile",'dns ss')
%%
%%
clf;
cd     'C:\Users\Artur\Documents\rotst\self-sim'

set(gcf,"Position",[   675   100   560   607])

r=0:0.1:10; reclw=(40./r).^2;

semilogy(r,reclw); ylim([10 3000]); set(gca,"YTick",[16,100:100:600,900:300:3000]); grid on; hold on;

R=[]; ress=[];
a=load("dataR10.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR8.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR7.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR6.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR5.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR4.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR3.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR2.5.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR2.mat"); R=[R;a.R]; ress=[ress;a.re0];
a=load("dataR1.7.mat"); R=[R;a.R]; ress=[ress;a.re0];




semilogy(R,ress,'-x');

set(gca,"XTick",sort([get(gca,"XTick"),1.6]))

