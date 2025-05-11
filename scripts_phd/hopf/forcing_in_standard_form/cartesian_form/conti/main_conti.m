mu=-0.015;
K=0.02;
beta=pi/4;

%%
KM=[];
muM=[];

for ii=1:1000
    main_forcing
    if (abs(r(end)-rex) < 1e-8)
        KM=[KM; K];
        muM=[muM; mu];
        mu=mu-1e-3

    else
        K=0.99*K
        K=1.01*K
        %     break;
    end
    ii
end
%%

KM=[];
muM=[];

for ii=1:1000
    main_forcing
    if (abs(r(end)-rex) <1e-8)
        KM=[KM; K];
        muM=[muM; mu];
        mu=mu-1e-3
        K=0.9*K
    else
        %         K=0.99*K
        K=1.01*K
        %     break;
    end
    ii
end
%%
KM1=KM;
muM1=muM;
%%
KM2=KM;
muM2=muM;
%%
KM3=KM;
muM3=muM;
%%
KM4=KM;
muM4=muM;
%%
KM5=KM;
muM5=muM;
%%
KM6=KM;
muM6=muM;
%%
plot(muM,KM)
%%
hold on;
plot(muM,KM)
plot(muM,(gamma*muM-1)/(cos(beta)+gamma*sin(beta)))
%%
hold on;
plot(muM1,KM1,'b-')
plot(muM2,KM2,'b-')
plot(muM3,KM3,'r-')
plot(muM4,KM4,'r-')
plot(muM5,KM5,'k-')
plot(muM6,KM6,'k-')
beta=pi/2;
plot(muM1,(gamma*muM1-1)/(cos(beta)+gamma*sin(beta))/2/pi,'b--')
beta=pi/4;
plot(muM1,(gamma*muM1-1)/(cos(beta)+gamma*sin(beta))/2/pi,'r--')
beta=pi/8*3;
plot(muM1,(gamma*muM1-1)/(cos(beta)+gamma*sin(beta))/2/pi,'k--')

xlabel("mu"); ylabel("K"); title("limiting K for stabilisation | subHopf gamma=-10");
grid on; grid minor; legend("beta = pi/2","beta = pi/2","beta = pi/4","beta = pi/4","beta = pi*3/8","beta = pi*3/8", ...
    "beta = pi/2 anal","beta = pi/4 anal","beta = pi*3/8 anal","Location","best")
%%
save("data.mat","KM6","KM5","KM4","KM3","KM2","KM1","muM6","muM5","muM4","muM3","muM2","muM1","-mat");