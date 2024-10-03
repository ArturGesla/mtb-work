%%
clf; 
set(gcf,"Position",[          24         267        1306         399]);
a1=load("resarrto148.mat"); a=a1.resarr(:,2:end);
a1=load("resarrto2000.mat"); a=[fliplr(a),a1.resarr(:,2:end)];

subplot(1,4,1:2); plot(a(1,:),a(4,:),'-'); title("Red vs Reh at abs inst");
xlabel("Reh"); ylabel("Red"); 
hold on;
plot(a(1,:),10*sqrt(a(1,:)),'-');
plot(a(1,:),a(1,:)./a(1,:)*24/sqrt(0.313),':k');
plot(1000,48.5,'r+');

legend("Red crit","Red edge of cavity","Bodewadt limit (LW97)","Serre et al. 2004","Location","best"); grid on;

axes('Position',[.17 .6 .05 .25]); 
plot(a(1,:),a(4,:),'-'); hold on; xlim([130 200]); ylim([0 150]);
plot(a(1,:),10*sqrt(a(1,:)),'-'); box on; grid on; 

axes('Position',[.25 .27 .05 .25]); 
plot(a(1,:),a(4,:),'-'); hold on; xlim([900 1100]); ylim([40 50]);
plot(a(1,:),10*sqrt(a(1,:)),'-'); box on; grid on; 
plot(a(1,:),a(1,:)./a(1,:)*24/sqrt(0.313),':k');
plot(1000,48.5,'r+');

subplot(1,4,3); plot(a(1,:),a(2,:)); title("ar vs Reh at abs inst");
xlabel("Reh"); ylabel("ar"); grid on;

subplot(1,4,4); plot(a(1,:),a(3,:)); title("ai vs Reh at abs inst");
xlabel("Reh"); ylabel("ai"); grid on;

