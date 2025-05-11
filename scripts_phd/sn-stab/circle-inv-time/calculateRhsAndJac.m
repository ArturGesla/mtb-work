function [g,jac]=calculateRhsAndJac(neq,nt,u,r)
%%
ii=[]; jj=[]; vv=[];

%u=[xr0,yr0,zr0,xr1,...,xi0,yi0,zi0,...]

%Lorenz
s=10;
b=8/3;
% r=24;
%   dxdt=s*(y-x)
%   dydt=x*(r-z)-y
%   dzdt=xy-b*z


lastid=nt*neq*2+1;
om=u(lastid);
 exchange=0;  for ik=0:2*nt-1
%      exchange=1;  for ik=[0,2:2*nt-1]
% exchange=1; for ik=[0:2:nt,nt:2:2*nt-1]%0:2*nt-1
% exchange=1; for ik=[0:nt-1,nt+1:2*nt-1]
    ix=ik*neq+1; iy=ix+1; iz=iy+1;
    
    realEq=ik<nt;



    kSum=mod(ik,nt);

    
         g(ix)=0;
         g(iy)=0;
         g(iz)=0;

    
%linear terms
for ikl=-nt+1:nt-1

    ixlr=abs(ikl)*neq+1; ixli=(abs(ikl)+nt)*neq+1; iylr=abs(ikl)*neq+2; iyli=(abs(ikl)+nt)*neq+2; izlr=abs(ikl)*neq+3; izli=(abs(ikl)+nt)*neq+3; %indices

if(realEq==1)
    g(ix)=g(ix)-1*(cr(ikl,kSum)*u(iylr)-ci(ikl,kSum,nt)*sgn(ikl)*u(iyli));
    g(ix)=g(ix)-ikl*om*(-u(ixlr)*ci_t(ikl,kSum,nt)-sgn(ikl)*u(ixli)*cr_t(ikl,kSum));
%     g(ix)

    g(iy)=g(iy)+1*(cr(ikl,kSum)*u(ixlr)-ci(ikl,kSum,nt)*sgn(ikl)*u(ixli));
    g(iy)=g(iy)-ikl*om*(-u(iylr)*ci_t(ikl,kSum,nt)-sgn(ikl)*u(iyli)*cr_t(ikl,kSum));

    g(iz)=g(iz)-1*(cr(ikl,kSum)*u(izlr)-ci(ikl,kSum,nt)*sgn(ikl)*u(izli));
    

%     %temp jac
%     ii(end+1)=ix; jj(end+1)=lastid; vv(end+1)=-ikl*(-u(ixlr)*ci_t(ikl,kSum,nt)-sgn(ikl)*u(ixli)*cr_t(ikl,kSum));
%     ii(end+1)=iy; jj(end+1)=lastid; vv(end+1)=-ikl*(-u(iylr)*ci_t(ikl,kSum,nt)-sgn(ikl)*u(iyli)*cr_t(ikl,kSum));
%     ii(end+1)=iz; jj(end+1)=lastid; vv(end+1)=-ikl*(-u(izlr)*ci_t(ikl,kSum,nt)-sgn(ikl)*u(izli)*cr_t(ikl,kSum));
%  
%     %jac
%     ii(end+1)=ix; jj(end+1)=ixlr; vv(end+1)=-s*(cr(ikl,kSum))               -ikl*om*(-ci_t(ikl,kSum,nt));
%     ii(end+1)=ix; jj(end+1)=ixli; vv(end+1)=-s*(-ci(ikl,kSum,nt)*sgn(ikl))    -ikl*om*(-sgn(ikl)*cr_t(ikl,kSum));
%     ii(end+1)=ix; jj(end+1)=iylr; vv(end+1)=s*(cr(ikl,kSum));
%     ii(end+1)=ix; jj(end+1)=iyli; vv(end+1)=s*(-ci(ikl,kSum,nt)*sgn(ikl));
% 
%     ii(end+1)=iy; jj(end+1)=iylr; vv(end+1)=-1*(cr(ikl,kSum))               -ikl*om*(-ci_t(ikl,kSum,nt));
%     ii(end+1)=iy; jj(end+1)=iyli; vv(end+1)=-1*(-ci(ikl,kSum,nt)*sgn(ikl))    -ikl*om*(-sgn(ikl)*cr_t(ikl,kSum));
%     ii(end+1)=iy; jj(end+1)=ixlr; vv(end+1)=r*(cr(ikl,kSum));
%     ii(end+1)=iy; jj(end+1)=ixli; vv(end+1)=r*(-ci(ikl,kSum,nt)*sgn(ikl));
% 
%     ii(end+1)=iz; jj(end+1)=izlr; vv(end+1)=-b*(cr(ikl,kSum))                 -ikl*om*(-ci_t(ikl,kSum,nt));
%     ii(end+1)=iz; jj(end+1)=izli; vv(end+1)=-b*(-ci(ikl,kSum,nt)*sgn(ikl))      -ikl*om*(-sgn(ikl)*cr_t(ikl,kSum));


else
    g(ix)=g(ix)-1*(ci(ikl,kSum,nt)*u(iylr)+cr(ikl,kSum)*sgn(ikl)*u(iyli));
    g(ix)=g(ix)-ikl*om*(+u(ixlr)*cr_t(ikl,kSum)-sgn(ikl)*u(ixli)*ci_t(ikl,kSum,nt));

    g(iy)=g(iy)+1*(ci(ikl,kSum,nt)*u(ixlr)+cr(ikl,kSum)*sgn(ikl)*u(ixli));
    g(iy)=g(iy)-ikl*om*(+u(iylr)*cr_t(ikl,kSum)-sgn(ikl)*u(iyli)*ci_t(ikl,kSum,nt));

    g(iz)=g(iz)-1*(ci(ikl,kSum,nt)*u(izlr)+cr(ikl,kSum)*sgn(ikl)*u(izli));

    %temp jac
%     ii(end+1)=ix; jj(end+1)=lastid; vv(end+1)=-ikl*(+u(ixlr)*cr_t(ikl,kSum)-sgn(ikl)*u(ixli)*ci_t(ikl,kSum,nt));
%     ii(end+1)=iy; jj(end+1)=lastid; vv(end+1)=-ikl*(+u(iylr)*cr_t(ikl,kSum)-sgn(ikl)*u(iyli)*ci_t(ikl,kSum,nt));
%     ii(end+1)=iz; jj(end+1)=lastid; vv(end+1)=-ikl*(+u(izlr)*cr_t(ikl,kSum)-sgn(ikl)*u(izli)*ci_t(ikl,kSum,nt));
%  
%     %jac
%     ii(end+1)=ix; jj(end+1)=ixlr; vv(end+1)=-s*(ci(ikl,kSum,nt))               -ikl*om*(cr_t(ikl,kSum));
%     ii(end+1)=ix; jj(end+1)=ixli; vv(end+1)=-s*(cr(ikl,kSum)*sgn(ikl))     -ikl*om*(-sgn(ikl)*ci_t(ikl,kSum,nt));
%     ii(end+1)=ix; jj(end+1)=iylr; vv(end+1)=s*(ci(ikl,kSum,nt));
%     ii(end+1)=ix; jj(end+1)=iyli; vv(end+1)=s*(cr(ikl,kSum)*sgn(ikl));
% 
%     ii(end+1)=iy; jj(end+1)=iylr; vv(end+1)=-1*(ci(ikl,kSum,nt))              -ikl*om*(cr_t(ikl,kSum));
%     ii(end+1)=iy; jj(end+1)=iyli; vv(end+1)=-1*(cr(ikl,kSum)*sgn(ikl))    -ikl*om*(-sgn(ikl)*ci_t(ikl,kSum,nt));
%     ii(end+1)=iy; jj(end+1)=ixlr; vv(end+1)=+r*(ci(ikl,kSum,nt));
%     ii(end+1)=iy; jj(end+1)=ixli; vv(end+1)=+r*(cr(ikl,kSum)*sgn(ikl));
% 
%     ii(end+1)=iz; jj(end+1)=izlr; vv(end+1)=-b*(ci(ikl,kSum,nt))                 -ikl*om*(cr_t(ikl,kSum));
%     ii(end+1)=iz; jj(end+1)=izli; vv(end+1)=-b*(cr(ikl,kSum)*sgn(ikl))       -ikl*om*(-sgn(ikl)*ci_t(ikl,kSum,nt));
end

end
    
end

if(false)%(exchange==1)

%     %then the half mode equations without first real mode
% for ik=[3:2:nt-1,nt+1:2:2*nt-1]
%     ix=ik*neq+1; iy=ix+1; iz=iy+1;
%          g(ix)=u(ix);
%          g(iy)=u(iy);
%          g(iz)=u(iz);
% ii(end+1)=ix; jj(end+1)=ix; vv(end+1)=1;
% ii(end+1)=iy; jj(end+1)=iy; vv(end+1)=1;
% ii(end+1)=iz; jj(end+1)=iz; vv(end+1)=1;
% end

%finally the real half mode equation
for ik=1%nt
% for ik=[1:2:nt-1,nt+1:2:2*nt-1]
    ix=ik*neq+1; iy=ix+1; iz=iy+1;
    
%     has to express that u(1)-u(0)=0
for ikl=-nt+1+1:2:nt-1-1 %half mode only
    ixlr=abs(ikl)*neq+1; ixli=(abs(ikl)+nt)*neq+1; iylr=abs(ikl)*neq+2; iyli=(abs(ikl)+nt)*neq+2; izlr=abs(ikl)*neq+3; izli=(abs(ikl)+nt)*neq+3; %indices
    
    
    g(ix)=g(ix)-u(ixlr)-u(ixlr);
    g(iy)=g(iy)-u(iylr)-u(iylr);
    g(iz)=g(iz)-u(izlr)-u(izlr);
    
    ii(end+1)=ix; jj(end+1)=ixlr; vv(end+1)=-2;
ii(end+1)=iy; jj(end+1)=iylr; vv(end+1)=-2;
ii(end+1)=iz; jj(end+1)=izlr; vv(end+1)=-2;

         
end
         
end
end

%last eq %phase in 0
% g(lastid)=sum(u(neq*nt+1:neq:end-1).*[0:1:nt-1]');
% aa=neq*nt+1:neq:lastid-1;
% ii=[ii,aa./aa*lastid]; jj=[jj,aa]; vv=[vv,[0:1:nt-1]];

%last eq %phase in phi
phi=   0.433132291618699;
g(lastid)=sum(u(1:neq:neq*nt).*[0:1:nt-1]'.*sin(phi.*[0:1:nt-1]')+u(neq*nt+1:neq:end-1).*[0:1:nt-1]'.*cos(phi.*[0:1:nt-1]')); %time der at phi
aar=1:neq:neq*nt;
aai=neq*nt+1:neq:lastid-1;
ii=[ii,aar./aar*lastid]; jj=[jj,aar]; vv=[vv,([0:1:nt-1]'.*sin(phi.*[0:1:nt-1]'))'];
ii=[ii,aai./aai*lastid]; jj=[jj,aai]; vv=[vv,([0:1:nt-1]'.*cos(phi.*[0:1:nt-1]'))'];


jac=sparse(ii,jj,vv);
end