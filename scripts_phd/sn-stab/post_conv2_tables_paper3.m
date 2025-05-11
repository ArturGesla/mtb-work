%% rho =24;
clc;

% for i=[12,14,22,32,52]
%     for i=[32,52,72,92,112]-2
        for i=[52,72,92,112,132,152]-2
% a=load("lorenz-cheb-coll\flnum-"+num2str(i)+"-24-cheb-coll-lornez.mat");  
% a=load("lorenz-cheb-coll\flnum-"+num2str(i)+"-28-cheb-coll-lornez.mat");  
a=load("lorenz-cheb-coll\flnum-"+num2str(i)+"-160-cheb-coll-lornez.mat");  
fl=sort(real(a.exponents)); fl1=floor(log10(abs([fl,2*pi/a.u(end)]))); fl2=[fl,2*pi/a.u(end)]./10.^fl1;
fprintf("$%d$",a.nt-1);
fprintf("& $%1.4f \\times 10^{%d}$",fl2(1),fl1(1)); 
fprintf("& $%1.4f \\times 10^{%d}$",fl2(2),fl1(2));
fprintf("& $%1.4f \\times 10^{%d}$",fl2(3),fl1(3)); 
if (fl1(4)==0)
    fprintf("& $%1.8f$",fl2(4)); 
else 
fprintf("& $%1.8f \\times 10^{%d}$",fl2(4),fl1(4)); 
end
fprintf("\\\\ \n");
end

