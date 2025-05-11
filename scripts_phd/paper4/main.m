addpath     'C:\Users\Artur\Documents\GitHub\rotst2\scripts\source_for_mtb'


% Rec 2.2 104 104 1.2 104 9000 8700 7740 8500 10500 1.7 104 5.9 104 6 104
% σi 0.0108 0.0422 0.216 0.227 0.239 0.227 0.220 0.218 0.218 0.219 0.392

%%
del =[0.99 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.0];
si=[0.0108 0.0422 0.216 0.227 0.239 0.227 0.220 0.218 0.218 0.219 0.392];
re=[2.2e4 1e4 1.2e4 9000 8700 7740 8500 10500 1.7e4 5.9e4 6e4];

%%
clf;
plot(del,re,'-x'); grid on; xlabel("$R_2/R_1$"); ylabel("$Re_c$");
fnts=12; jfm_plt_aid_comm;


axes('position', [0.52 0.52 0.3 0.3]); 
plot(del,si,'-x'); xlabel("$R_2/R_1$"); ylabel("$\sigma_i$"); grid on; grid minor;
fnts=12; jfm_plt_aid_comm;

%%

