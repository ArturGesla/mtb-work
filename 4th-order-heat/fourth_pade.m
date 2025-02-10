A=zeros(length(xc)*2);
b=zeros(length(xc)*2,1);


%for each internal node
for i=2:length(xc)-1

    %pade compact
    iu=(i-1)*2+1;
    idu=(i-1)*2+1+1;
    iuxp=(i-1)*2+1+2;
    iduxp=(i-1)*2+1+3;
    iuxm=(i-1)*2+1-2;
    iduxm=(i-1)*2+1+1-2;

    A(iu,iduxp)=1/10;
    A(iu,idu)=1;
    A(iu,iduxm)=1/10;
    A(iu,iuxp)=-6/5/dx/dx;
    A(iu,iu)=+12/5/dx/dx;
    A(iu,iuxm)=-6/5/dx/dx;
    b(iu)=0;

    %constarin
    A(idu,idu)=1;
    b(idu)=f(xc(i));



end
%%
A(1,1)=[1]; A(end-1,end-1)=[1 ]; % hom bc

    A(2,2)=1; b(2)=f(xc(1));
    A(end,end)=1; b(end)=f(xc(end));
