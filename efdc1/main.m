% cd     '/home/gesla/Downloads'
clear;
a=importdata("Feuille de calcul sans titre - Arkusz1(1).csv");
addpath     '/home/gesla/Documents/git/rotst2/scripts/source_for_mtb';

%%
for i=1:length(a.data)
num(i)=a.data(i,1);
ysa(i)=a.data(i,2);
% if(min(abs(num(i)-num(1:i-1)))==0)
%     num(i)=num(i)+1;
%     i=i-1
% end
lbl(i)=(a.textdata(i+1,1));
multin(i)=mod(str2double(a.textdata(i+1,2)),10);
re(i)=str2double(a.textdata(i+1,3));
desc(i)=(a.textdata(i+1,4));
desc(i)={strrep(desc(i)," ","\n")};
end
%%

% clf;
% for i=1:length(num)
% bubblechart(multin(i),re(i),sqrt(num(i))); hold on;
% b=text(multin(i),re(i),string(lbl(i))+desc(i)+" ("+num2str(num(i))+")",'HorizontalAlignment','center','FontSize',sqrt(num(i))*1.5);
% end
% fnts=12; jfm_plt_aid_comm;
% bubblesize([min(sqrt(num)) max(sqrt(num))]*10)
% %%
% clf;
% 
% bubblechart([1 2 ],[2 2],[1 2]);
% bubblesize([1 2]*100); grid on; grid minor;

%%
clc;
for i=1:length(num)
    fprintf("%d %f\n",1,sqrt(num(i)));
end
%%
%http://hydra.nat.uni-magdeburg.de/packing/wizard/wizard.php

a=importdata('wizard-Zu9YMG-34.coords');
% a1=abs(a.data(:,3)-round(sqrt(num),6));
% [a2,b]=min(a1);

% [a1,b1]=sort(a.data(:,3));
[a2,b2]=sort(sqrt(num));


num2=num(b2);
lbl2=lbl(b2);
desc2=desc(b2);
ysa2=ysa(b2);
%%
for i=1:34
desc3(i)=sprintf(string(desc2(i)));
end
%%


%
clarr=["4285f4","db4437","f4b400","0f9d58","ab47bc","00acc1","ff7043","9e9d24","5c6bc0","f06292","00796b","ca3771","b3cefb","f7b4ae","fde49b","aedcba","ffc599","b5e5e8","ecf3fe","fdeceb","fff8e6","ebf6ee","fff9f6","edf8f9"];
%
% desc2{34}=


clf;

 axis equal;
for i=1:length(a.data)
% bubblechart(a.data(i,1),a.data(i,2),a.data(i,3)); hold on;

str="#"+clarr(1+mod(i-1,18));
circle(a.data(i,1),a.data(i,2),a.data(i,3)/2,str,0.6); hold on;

% str="#ff0000";
% circle(a.data(i,1),a.data(i,2),a.data(i,3)/2,str,ysa2(i)/num2(i)); hold on;


% b=text(a.data(i,1),a.data(i,2),string(lbl2(i))+desc2(i)+" ("+num2str(num2(i))+")",'HorizontalAlignment','center','FontSize',sqrt(num2(i))*1.5);
b=text(a.data(i,1),a.data(i,2),desc3(i)+" ("+num2str(num2(i))+")",'HorizontalAlignment','center','FontSize',a.data(i,3)*1.7);
end
% fnts=12; jfm_plt_aid_comm; set(gcf,"Position",[   515   369   570   413])
% bubblesize([min(a.data(:,3)) max(a.data(:,3))]*10.2);
axis off;
title("EFDC1 Aachen (total submissions: "+num2str(sum(num))+")")

%%

clf; hold on; axis equal;  
circle(1,1,1)
circle(1,3,2)