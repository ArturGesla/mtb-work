function [g,jac]=calculateRhsAndJac(neq,nt,u,lam,om1,collx)
%%

%s=10;
%b=8/3;
% r=24;
% lam=r.^2;
ii=[]; jj=[]; vv=[];

g=u*0;

%u=[xr0,yr0,zr0,xr1,...,xi0,yi0,zi0,...]

% Noack system
% f=[mu*u-v*g-w*u;
%     mu*v+u/g-v*w;
%     -w+u^2+v^2*g*g];

lastid=nt*neq+1;
om=u(lastid);

%non-linear
for ikl=0:nt-1
    for ikr=0:nt-1

        ixr=abs(ikr)*neq+1;  ixl=abs(ikl)*neq+1; iyr=abs(ikr)*neq+2;  iyl=abs(ikl)*neq+2; izr=abs(ikr)*neq+3;  izl=abs(ikl)*neq+3;

        g(1:neq:end-1)=g(1:neq:end-1)+1/2*(u(ixl)*u(izr))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)));
        g(2:neq:end-1)=g(2:neq:end-1)+1/2*(u(izl)*u(iyr))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)));
        g(3:neq:end-1)=g(3:neq:end-1)-1/2*(u(ixl)*u(ixr))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)))-1/2*(u(iyl)*u(iyr))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)))-1/2*(u(izl)*u(izr))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)));

        ii=[ii,1:neq:nt*neq]; jj=[jj,ixl*ones([1,nt])];  vv=[vv; -1/2*(u(izr))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)))];
        ii=[ii,1:neq:nt*neq]; jj=[jj,izr*ones([1,nt])];  vv=[vv; -1/2*(u(ixl))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)))];

        ii=[ii,2:neq:nt*neq]; jj=[jj,izl*ones([1,nt])];  vv=[vv; -1/2*(u(iyr))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)))];
        ii=[ii,2:neq:nt*neq]; jj=[jj,iyr*ones([1,nt])];  vv=[vv; -1/2*(u(izl))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)))];

        ii=[ii,3:neq:nt*neq]; jj=[jj,ixl*ones([1,nt])];  vv=[vv; 1/2*(u(ixr))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)))];
        ii=[ii,3:neq:nt*neq]; jj=[jj,ixr*ones([1,nt])];  vv=[vv; 1/2*(u(ixl))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)))];
        ii=[ii,3:neq:nt*neq]; jj=[jj,iyr*ones([1,nt])];  vv=[vv; 1/2*(u(iyl))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)))];
        ii=[ii,3:neq:nt*neq]; jj=[jj,iyl*ones([1,nt])];  vv=[vv; 1/2*(u(iyr))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)))];

    end
end

%triple nonlinearity
for ikl=0:nt-1
    for ikr=0:nt-1
        for ikm=0:nt-1

        ixr=abs(ikr)*neq+1;  ixl=abs(ikl)*neq+1; iyr=abs(ikr)*neq+2;  iyl=abs(ikl)*neq+2; izr=abs(ikr)*neq+3;  izl=abs(ikl)*neq+3; ixm=abs(ikm)*neq+1;   iym=abs(ikm)*neq+2;  izm=abs(ikm)*neq+3; 

        %no projection of product to cheb modes - actually not so usefull on collocation
        g(1:neq:end-1)=g(1:neq:end-1)-(u(ixl)*u(izr)*u(izm))*(cos((ikl)*acos(collx)).*cos((ikr)*acos(collx)).*cos((ikm)*acos(collx)));
        g(2:neq:end-1)=g(2:neq:end-1)-(u(iyl)*u(izr)*u(izm))*(cos((ikl)*acos(collx)).*cos((ikr)*acos(collx)).*cos((ikm)*acos(collx)));
        % g(3:neq:end-1)=g(3:neq:end-1)-1/2*(u(ixl)*u(ixr))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)))-1/2*(u(iyl)*u(iyr))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)))-1/2*(u(izl)*u(izr))*(cos((ikl+ikr)*acos(collx))+cos(abs(ikl-ikr)*acos(collx)));
        end
    end
end



%linear part
for ikl=0:nt-1

    ixr=abs(ikr)*neq+1;  ixl=abs(ikl)*neq+1; iyr=abs(ikr)*neq+2;  iyl=abs(ikl)*neq+2; izr=abs(ikr)*neq+3;  izl=abs(ikl)*neq+3;

    g(1:neq:end-1)=g(1:neq:end-1)+(0.2+lam-3)*(u(ixl))*(cos((ikl)*acos(collx)))-0.25*(u(iyl))*(cos((ikl)*acos(collx)));
    g(2:neq:end-1)=g(2:neq:end-1)+(lam-3+0.2)*(u(iyl))*(cos((ikl)*acos(collx)))+0.25*(u(ixl))*(cos((ikl)*acos(collx)));
    g(3:neq:end-1)=g(3:neq:end-1)+lam*(u(izl))*(cos((ikl)*acos(collx)));


    ii=[ii,1:neq:nt*neq]; jj=[jj,ixl*ones([1,nt])];  vv=[vv; (0.2+lam-3)*(cos((ikl)*acos(collx)))];
    ii=[ii,1:neq:nt*neq]; jj=[jj,iyl*ones([1,nt])];  vv=[vv; -0.25*(cos((ikl)*acos(collx)))];

    ii=[ii,2:neq:nt*neq]; jj=[jj,iyl*ones([1,nt])];  vv=[vv; (0.2+lam-3)*(cos((ikl)*acos(collx)))];
    ii=[ii,2:neq:nt*neq]; jj=[jj,ixl*ones([1,nt])];  vv=[vv; +0.25*(cos((ikl)*acos(collx)))];

    ii=[ii,3:neq:nt*neq]; jj=[jj,izl*ones([1,nt])];  vv=[vv; lam*(cos((ikl)*acos(collx)))];

end
%
% % temporal term
for ikr=0:nt-1
    for ikl=mod(ikr+1,2):2:ikr

        ixr=abs(ikr)*neq+1; iyr=abs(ikr)*neq+2; izr=abs(ikr)*neq+3;

        g(1:neq:end-1)=g(1:neq:end-1)-ikr*(2*cos((ikl)*acos(collx))-1*(ikl==0))*om/pi/om1*u(ixr);
        g(2:neq:end-1)=g(2:neq:end-1)-ikr*(2*cos((ikl)*acos(collx))-1*(ikl==0))*om/pi/om1*u(iyr);
        g(3:neq:end-1)=g(3:neq:end-1)-ikr*(2*cos((ikl)*acos(collx))-1*(ikl==0))*om/pi/om1*u(izr);

        ii=[ii,1:neq:nt*neq]; jj=[jj,ixr*ones([1,nt])];  vv=[vv; -ikr*(2*cos((ikl)*acos(collx))-1*(ikl==0))*om/pi/om1];
        ii=[ii,1:neq:nt*neq]; jj=[jj,lastid*ones([1,nt])];  vv=[vv; -ikr*(2*cos((ikl)*acos(collx))-1*(ikl==0))/pi/om1*u(ixr)];

        ii=[ii,2:neq:nt*neq]; jj=[jj,iyr*ones([1,nt])];  vv=[vv; -ikr*(2*cos((ikl)*acos(collx))-1*(ikl==0))*om/pi/om1];
        ii=[ii,2:neq:nt*neq]; jj=[jj,lastid*ones([1,nt])];  vv=[vv; -ikr*(2*cos((ikl)*acos(collx))-1*(ikl==0))/pi/om1*u(iyr)];

        ii=[ii,3:neq:nt*neq]; jj=[jj,izr*ones([1,nt])];  vv=[vv; -ikr*(2*cos((ikl)*acos(collx))-1*(ikl==0))*om/pi/om1];
        ii=[ii,3:neq:nt*neq]; jj=[jj,lastid*ones([1,nt])];  vv=[vv; -ikr*(2*cos((ikl)*acos(collx))-1*(ikl==0))/pi/om1*u(izr)];




    end
end

%     for ikl=ik+1:2:nt-1
%     ixl=abs(ikl)*neq+1; iyl=abs(ikl)*neq+2; izl=abs(ikl)*neq+3; 
%   
%     g(ix)=g(ix)-u(ixl)*ikl*(2-(1*(ik==0)))*om/pi*1/om1;
%     g(iy)=g(iy)-u(iyl)*ikl*(2-(1*(ik==0)))*om/pi*1/om1;
%     g(iz)=g(iz)-u(izl)*ikl*(2-(1*(ik==0)))*om/pi;
%     
%     ii(end+1)=ix; jj(end+1)=ixl; vv(end+1)=-ikl*(2-(1*(ik==0)))*om/pi*1/om1;
%     ii(end+1)=iy; jj(end+1)=iyl; vv(end+1)=-ikl*(2-(1*(ik==0)))*om/pi*1/om1;
%     ii(end+1)=iz; jj(end+1)=izl; vv(end+1)=-ikl*(2-(1*(ik==0)))*om/pi;
%     
%     ii(end+1)=ix; jj(end+1)=lastid; vv(end+1)=-u(ixl)*ikl*(2-(1*(ik==0)))/pi*1/om1;
%     ii(end+1)=iy; jj(end+1)=lastid; vv(end+1)=-u(iyl)*ikl*(2-(1*(ik==0)))/pi*1/om1;
%     ii(end+1)=iz; jj(end+1)=lastid; vv(end+1)=-u(izl)*ikl*(2-(1*(ik==0)))/pi;
%                 
%     end
%     
%     
%     
% end
% g(lastid)=sum(u(neq*nt+1:neq:end-1).*[0:1:nt-1]');
% aa=neq*nt+1:neq:lastid-1;
% ii=[ii,aa./aa*lastid]; jj=[jj,aa]; vv=[vv,[0:1:nt-1]];

%evp
% for ik=nt-1
%     
%     ix=ik*neq+1; iy=ix+1; iz=iy+1;
%     
%     %has to express u(1)-u(0)=0
%     for ikl=0:nt-1
%     ixl=abs(ikl)*neq+1; iyl=abs(ikl)*neq+2; izl=abs(ikl)*neq+3; 
%   
%     g(ix)=g(ix)+u(ixl)-(-1)^ikl*u(ixl);
%     g(iy)=g(iy)+u(iyl)-(-1)^ikl*u(iyl);
%     g(iz)=g(iz)+u(izl)-(-1)^ikl*u(izl);
%     
%     ii(end+1)=ix; jj(end+1)=ixl; vv(end+1)=(1-(-1)^ikl);
%     ii(end+1)=iy; jj(end+1)=iyl; vv(end+1)=(1-(-1)^ikl);
%     ii(end+1)=iz; jj(end+1)=izl; vv(end+1)=(1-(-1)^ikl);
%     
%     end
% end
    






    
%phase
for ikl=1:nt-1
    ixl=abs(ikl)*neq+1; %iyl=abs(ikl)*neq+2; izl=abs(ikl)*neq+3; 
  
    g(lastid)=g(lastid)+u(ixl)*ikl^2; 
    
    ii(end+1)=lastid; jj(end+1)=ixl; vv(end+1)=ikl^2; 
                
end

ii(end+1)=lastid;jj(end+1)=lastid;vv(end+1)=0;
jac=sparse(ii,jj,vv);

%exchange
jac(end-3:end-1,:)=0; g(end-3:end-1)=0;

for ik=nt-1
    
    ix=ik*neq+1; iy=ix+1; iz=iy+1;
    
    %has to express u(1)-u(0)=0
    for ikl=0:nt-1
    ixl=abs(ikl)*neq+1; iyl=abs(ikl)*neq+2; izl=abs(ikl)*neq+3; 
  
    g(ix)=g(ix)+u(ixl)-(-1)^ikl*u(ixl);
    g(iy)=g(iy)+u(iyl)-(-1)^ikl*u(iyl);
    g(iz)=g(iz)+u(izl)-(-1)^ikl*u(izl);
    
%     ii(end+1)=ix; jj(end+1)=ixl; vv(end+1)=(1-(-1)^ikl);
%     ii(end+1)=iy; jj(end+1)=iyl; vv(end+1)=(1-(-1)^ikl);
%     ii(end+1)=iz; jj(end+1)=izl; vv(end+1)=(1-(-1)^ikl);

jac(ix,ixl)=jac(ix,ixl)+(1-(-1)^ikl);
jac(iy,iyl)=jac(iy,iyl)+(1-(-1)^ikl);
jac(iz,izl)=jac(iz,izl)+(1-(-1)^ikl);
    
    end
end


end