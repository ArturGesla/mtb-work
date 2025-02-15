clf;
set(gcf,"Position",[                46         381        1113         401])

tiledlayout(1,2); 
nexttile;

x=-cos(linspace(0,pi,2*nx))'; y=x;
% x=(linspace(-1,1,nx))'; y=x;
% u=real(ev(:,iev)); iev=iev+1;
% uarr=[uarr;u((length(x)+1)/2)]
uPhys=zeros(length(x)*length(y),1);
for ix=1:length(x)
    for iy=1:length(y)
        ip=iy+(ix-1)*length(y); xc=x(ix); yc=y(iy);

        for ikx=0:nt-1
            for iky=0:nt-1
                ipk=iky+1+(ikx-1+1)*ny;
                uPhys(ip)=uPhys(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
            end
        end
    end
end
%
uPhys2=reshape(uPhys,[length(y),length(x)]);
mesh(x,y,uPhys2);
xlabel('x'); ylabel('y');
title("\Delta u=0")

nexttile; 

% semilogy(narr,abs(uarr-2/pi/pi),'x-'); hold on;
% semilogy(narr,abs(uarr-uarr(end)),'x-'); hold on;
loglog(narr,abs(uarr-uarr(end,:)),'x-'); hold on;
% semilogy(narr,narr.^(-2),'o-')
% loglog(narr,narr.^(-2),'o-')
% loglog(narr,narr.^(-3),'o-')
% loglog(narr,narr.^(-4),'o-')
% loglog(narr,narr.^(-5),'o-')
% loglog(narr,narr.^(-12),'o-')

%
err=abs(uarr-uarr(end,:));

a=polyfit(log(narr(1:end-4)),log(err(1:end-4)),1);
loglog(narr,exp(a(1)*log(narr)+a(2)),'-'); 
legend("error |un(0.1,0.1)-ulast(0.1,0.1)|","fit y="+num2str(a(1)+"x")); grid on;
xlabel('n'); ylabel('error');
title("error of Cheb solution")
%%
exportgraphics(gcf,"chebStudy-rotst-0.03.png")

