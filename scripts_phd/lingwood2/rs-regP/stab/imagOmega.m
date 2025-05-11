function imOmg=imagOmega(ar,ai,bbar,R,data)

% a=load("../vk-np-130-k-1.mat");
% a=load("../rs-np-132-k-0.313-L-32.mat");

a=data;

x=a.x;
u=a.u*0;
U=a.u;


% R=530; beta=67; bbar=beta/R; 
omega=0;


% alpha=ar(ir)+1i*ai(ii);
alpha=ar+1i*ai;
     [g,jac0,jac1,jac2,jacom]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);
    jac=jac0+alpha*jac1+alpha.^2*jac2;
    [ev]=eigs(jac,-jacom,5,-0.2+0.1i);
%     eva=[eva,ev];
%     disp(alpha);
% fprintf("ir %d\t ii %d\n",ir,ii);
    [a,b]=max(imag(ev));
%     z(ir,ii)=ev(b);
    imOmg=(ev(b));
    
end