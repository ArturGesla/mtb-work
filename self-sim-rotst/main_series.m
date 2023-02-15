sqrt(u(end))
%%
al=0.34;
rl=21.6;
evsM={};

for al=0:0.0025:0.5
    evsMM=[];
    for rl=[17 18.2658 21.6]
%         for rl=2%:4:22
%             for rl=21.6%18.44%:4:22
%         k=0.313;
        k=sqrt(u(end));
        delta=L/sqrt(Re*k);
        kr=al/delta;
        ra=rl*delta;

        [rhs,jac,B]=calculateJacAndRhs(zc,zw,u,Re,ra,kr,L);
%         !calculateJacAndRhs "zc" "zw" "u" "Re" "ra" "kr" "L"
        [evc,evs]=eigs((jac),(B),40,"smallestabs"); evs=diag(evs);
        evsMM=[evsMM,evs];
        fprintf("al: %4.2e \t rl: %4.2f \n",al,rl);
    end
    evsM{end+1}=evsMM;
end

%%
al=0.2775;
rl=18.2658;
k=sqrt(u(end));
        delta=L/sqrt(Re*k);
        kr=al/delta;
        ra=rl*delta;
[rhs,jac,B]=calculateJacAndRhs(zc,zw,u,Re,ra,kr,L);
        [evc,evs]=eigs((jac),(B),80,"smallestabs"); evs=diag(evs);

[c,b]=max(real(evs))
%%
close all;
hold on; grid on; grid minor;
symbol=repmat(".",1,length(evsM));
% symbol(1)="x";
symbol(112)="o";
colour=["b","r","g"];
for i=1:length(evsM)
plot(evsM{i},symbol(i)); set(gca,'ColorOrderIndex',1);
end
xlim([-1 0.05]);
xlabel("real"); ylabel("imag");
title("Re: "+num2str(Re));
exportgraphics(gcf,"spectrum2-"+num2str(Re)+".png","Resolution",150);
%%
alM=0:0.0025:0.5;
alM(112)
%%
ire=2; maxlam=-10;
for i=1:length(evsM)
if (max(real((evsM{i}(:,ire))))>maxlam)
    [maxlam,f]=max(real((evsM{i}(:,ire))));
%     f
    i
end

end
maxlam
%%
 interp1([-0.0029,0.0019],[18,18.44],[0])