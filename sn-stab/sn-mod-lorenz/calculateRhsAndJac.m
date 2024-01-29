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
for ik=0:2*nt-1
    ix=ik*neq+1; iy=ix+1; iz=iy+1;
    
    realEq=ik<nt;
    
%     
    nlxz=0;
    nlxy=0;
% 
%     nlxz=0;
%     nlyz=0;
%     nlyy=0;
%     nlxx=0;



    kSum=mod(ik,nt);

    %nonlinear term
%     for ikl=kSum-nt+1:nt-1
%         ikr=kSum-ikl;
    for ikl=-nt+1:nt-1
        for ikr=-nt+1:nt-1

        
        ixrr=abs(ikr)*neq+1; ixri=(abs(ikr)+nt)*neq+1; ixlr=abs(ikl)*neq+1; ixli=(abs(ikl)+nt)*neq+1;
        iyrr=abs(ikr)*neq+2; iyri=(abs(ikr)+nt)*neq+2; iylr=abs(ikl)*neq+2; iyli=(abs(ikl)+nt)*neq+2;
        izrr=abs(ikr)*neq+3; izri=(abs(ikr)+nt)*neq+3; izlr=abs(ikl)*neq+3; izli=(abs(ikl)+nt)*neq+3;
        
        if(realEq==1)
% ikl
% ikr
            nlxz=nlxz+cr(ikl+ikr,kSum)*(u(ixrr)*u(izlr)- sgn(ikr)* sgn(ikl)*u(ixli)*u(izri))-ci(ikl+ikr,kSum)*(u(ixrr)*u(izli)* sgn(ikl)+ sgn(ikr)*u(ixri)*u(izlr));
            nlxy=nlxy+cr(ikl+ikr,kSum)*(u(ixrr)*u(iylr)- sgn(ikr)* sgn(ikl)*u(ixli)*u(iyri))-ci(ikl+ikr,kSum)*(u(ixrr)*u(iyli)* sgn(ikl)+ sgn(ikr)*u(ixri)*u(iylr));
%             nlxx=nlxx+cr(ikl+ikr,kSum)*(u(ixrr)*u(ixlr)- sgn(ikr)* sgn(ikl)*u(ixli)*u(ixri))-ci(ikl+ikr,kSum)*(u(ixrr)*u(ixli)* sgn(ikl)+ sgn(ikr)*u(ixri)*u(ixlr));
%             nlyy=nlyy+cr(ikl+ikr,kSum)*(u(iyrr)*u(iylr)- sgn(ikr)* sgn(ikl)*u(iyli)*u(iyri))-ci(ikl+ikr,kSum)*(u(iyrr)*u(iyli)* sgn(ikl)+ sgn(ikr)*u(iyri)*u(iylr));

            %jac
            ii(end+1)=iy; jj(end+1)=ixrr; vv(end+1)=(+cr(ikl+ikr,kSum)*(u(izlr))-ci(ikl+ikr,kSum)*(u(izli)* sgn(ikl)))*-1;
            ii(end+1)=iy; jj(end+1)=izlr; vv(end+1)=(+cr(ikl+ikr,kSum)*(u(ixrr))-ci(ikl+ikr,kSum)*(sgn(ikr)*u(ixri)))*-1;
            ii(end+1)=iy; jj(end+1)=ixli; vv(end+1)=(+cr(ikl+ikr,kSum)*(- sgn(ikr)* sgn(ikl)*u(izri)))*-1;
            ii(end+1)=iy; jj(end+1)=izri; vv(end+1)=(+cr(ikl+ikr,kSum)*(- sgn(ikr)* sgn(ikl)*u(ixli)))*-1;
            ii(end+1)=iy; jj(end+1)=izli; vv(end+1)=(-ci(ikl+ikr,kSum)*(u(ixrr)* sgn(ikl)))*-1;
            ii(end+1)=iy; jj(end+1)=ixri; vv(end+1)=(-ci(ikl+ikr,kSum)*( sgn(ikr)*u(izlr)))*-1;
             
            ii(end+1)=iz; jj(end+1)=ixrr; vv(end+1)=(+cr(ikl+ikr,kSum)*(u(iylr))-ci(ikl+ikr,kSum)*(u(iyli)* sgn(ikl)))*1;
            ii(end+1)=iz; jj(end+1)=iylr; vv(end+1)=(+cr(ikl+ikr,kSum)*(u(ixrr))-ci(ikl+ikr,kSum)*( sgn(ikr)*u(ixri)))*1;
            ii(end+1)=iz; jj(end+1)=ixli; vv(end+1)=(+cr(ikl+ikr,kSum)*(- sgn(ikr)* sgn(ikl)*u(iyri)))*1;
            ii(end+1)=iz; jj(end+1)=iyri; vv(end+1)=(+cr(ikl+ikr,kSum)*(- sgn(ikr)* sgn(ikl)*u(ixli)))*1;
            ii(end+1)=iz; jj(end+1)=iyli; vv(end+1)=(-ci(ikl+ikr,kSum)*(u(ixrr)* sgn(ikl)))*1;
            ii(end+1)=iz; jj(end+1)=ixri; vv(end+1)=(-ci(ikl+ikr,kSum)*(sgn(ikr)*u(iylr)))*1;


            
            
        else

             nlxz=nlxz+ci(ikl+ikr,kSum)*(u(ixrr)*u(izlr)- sgn(ikr)* sgn(ikl)*u(ixli)*u(izri))+cr(ikl+ikr,kSum)*(u(ixrr)*u(izli)* sgn(ikl)+ sgn(ikr)*u(ixri)*u(izlr));
             nlxy=nlxy+ci(ikl+ikr,kSum)*(u(ixrr)*u(iylr)- sgn(ikr)* sgn(ikl)*u(ixli)*u(iyri))+cr(ikl+ikr,kSum)*(u(ixrr)*u(iyli)* sgn(ikl)+ sgn(ikr)*u(ixri)*u(iylr));
%              nlxx=nlxx+ci(ikl+ikr,kSum)*(u(ixrr)*u(ixlr)- sgn(ikr)* sgn(ikl)*u(ixli)*u(ixri))+cr(ikl+ikr,kSum)*(u(ixrr)*u(ixli)* sgn(ikl)+ sgn(ikr)*u(ixri)*u(ixlr));
%              nlyy=nlyy+ci(ikl+ikr,kSum)*(u(iyrr)*u(iylr)- sgn(ikr)* sgn(ikl)*u(iyli)*u(iyri))+cr(ikl+ikr,kSum)*(u(iyrr)*u(iyli)* sgn(ikl)+ sgn(ikr)*u(iyri)*u(iylr));

            %jac
            ii(end+1)=iy; jj(end+1)=ixrr; vv(end+1)=(+ci(ikl+ikr,kSum)*(u(izlr))+cr(ikl+ikr,kSum)*(u(izli)* sgn(ikl)))*-1;
            ii(end+1)=iy; jj(end+1)=izlr; vv(end+1)=(+ci(ikl+ikr,kSum)*(u(ixrr))+cr(ikl+ikr,kSum)*(sgn(ikr)*u(ixri)))*-1;
            ii(end+1)=iy; jj(end+1)=ixli; vv(end+1)=(+ci(ikl+ikr,kSum)*(- sgn(ikr)* sgn(ikl)*u(izri)))*-1;
            ii(end+1)=iy; jj(end+1)=izri; vv(end+1)=(+ci(ikl+ikr,kSum)*(- sgn(ikr)* sgn(ikl)*u(ixli)))*-1;
            ii(end+1)=iy; jj(end+1)=izli; vv(end+1)=(+cr(ikl+ikr,kSum)*(u(ixrr)* sgn(ikl)))*-1;
            ii(end+1)=iy; jj(end+1)=ixri; vv(end+1)=(+cr(ikl+ikr,kSum)*( sgn(ikr)*u(izlr)))*-1;

            ii(end+1)=iz; jj(end+1)=ixrr; vv(end+1)=(+ci(ikl+ikr,kSum)*(u(iylr))+cr(ikl+ikr,kSum)*(u(iyli)* sgn(ikl)))*1;
            ii(end+1)=iz; jj(end+1)=iylr; vv(end+1)=(+ci(ikl+ikr,kSum)*(u(ixrr))+cr(ikl+ikr,kSum)*( sgn(ikr)*u(ixri)))*1;
            ii(end+1)=iz; jj(end+1)=ixli; vv(end+1)=(+ci(ikl+ikr,kSum)*(- sgn(ikr)* sgn(ikl)*u(iyri)))*1;
            ii(end+1)=iz; jj(end+1)=iyri; vv(end+1)=(+ci(ikl+ikr,kSum)*(- sgn(ikr)* sgn(ikl)*u(ixli)))*1;
            ii(end+1)=iz; jj(end+1)=iyli; vv(end+1)=(+cr(ikl+ikr,kSum)*(u(ixrr)* sgn(ikl)))*1;
            ii(end+1)=iz; jj(end+1)=ixri; vv(end+1)=(+cr(ikl+ikr,kSum)*(sgn(ikr)*u(iylr)))*1;


%             ii(end+1)=iz; jj(end+1)=ixrr; vv(end+1)=+ci(ikl+ikr,kSum)*(u(ixlr))+cr(ikl+ikr,kSum)*(u(ixli)* sgn(ikl));
%             ii(end+1)=iz; jj(end+1)=ixlr; vv(end+1)=+ci(ikl+ikr,kSum)*(u(ixrr))+cr(ikl+ikr,kSum)*( sgn(ikr)*u(ixri));
%             ii(end+1)=iz; jj(end+1)=ixli; vv(end+1)=+ci(ikl+ikr,kSum)*(- sgn(ikr)* sgn(ikl)*u(ixri))+cr(ikl+ikr,kSum)*(u(ixrr)* sgn(ikl));
%             ii(end+1)=iz; jj(end+1)=ixri; vv(end+1)=+ci(ikl+ikr,kSum)*(- sgn(ikr)* sgn(ikl)*u(ixli))+cr(ikl+ikr,kSum)*(sgn(ikr)*u(ixlr));
% 
%             ii(end+1)=iz; jj(end+1)=iyrr; vv(end+1)=(ci(ikl+ikr,kSum)*(u(iylr))+cr(ikl+ikr,kSum)*(u(iyli)* sgn(ikl)))*gm*gm;
%             ii(end+1)=iz; jj(end+1)=iylr; vv(end+1)=(ci(ikl+ikr,kSum)*(u(iyrr))+cr(ikl+ikr,kSum)*( sgn(ikr)*u(iyri)))*gm*gm;
%             ii(end+1)=iz; jj(end+1)=iyli; vv(end+1)=(ci(ikl+ikr,kSum)*(- sgn(ikr)* sgn(ikl)*u(iyri))+cr(ikl+ikr,kSum)*(u(iyrr)* sgn(ikl)))*gm*gm;
%             ii(end+1)=iz; jj(end+1)=iyri; vv(end+1)=(ci(ikl+ikr,kSum)*(- sgn(ikr)* sgn(ikl)*u(iyli))+cr(ikl+ikr,kSum)*( sgn(ikr)*u(iylr)))*gm*gm;            
        end
        end
    end
    
         g(ix)=0;
         g(iy)=-nlxz;
         g(iz)=nlxy;

    
%linear terms
for ikl=-nt+1:nt-1

    ixlr=abs(ikl)*neq+1; ixli=(abs(ikl)+nt)*neq+1; iylr=abs(ikl)*neq+2; iyli=(abs(ikl)+nt)*neq+2; izlr=abs(ikl)*neq+3; izli=(abs(ikl)+nt)*neq+3; %indices

if(realEq==1)
    g(ix)=g(ix)+s*(cr(ikl,kSum)*u(iylr)-ci(ikl,kSum)*sgn(ikl)*u(iyli));
    g(ix)=g(ix)-s*(cr(ikl,kSum)*u(ixlr)-ci(ikl,kSum)*sgn(ikl)*u(ixli));
    g(ix)=g(ix)-ikl*om*(-u(ixlr)*ci(ikl,kSum)-sgn(ikl)*u(ixli)*cr(ikl,kSum));
%     g(ix)

    g(iy)=g(iy)+r*(cr(ikl,kSum)*u(ixlr)-ci(ikl,kSum)*sgn(ikl)*u(ixli));
    g(iy)=g(iy)-1*(cr(ikl,kSum)*u(iylr)-ci(ikl,kSum)*sgn(ikl)*u(iyli));
    g(iy)=g(iy)-ikl*om*(-u(iylr)*ci(ikl,kSum)-sgn(ikl)*u(iyli)*cr(ikl,kSum));

    g(iz)=g(iz)-b*(cr(ikl,kSum)*u(izlr)-ci(ikl,kSum)*sgn(ikl)*u(izli));
    g(iz)=g(iz)-ikl*om*(-u(izlr)*ci(ikl,kSum)-sgn(ikl)*u(izli)*cr(ikl,kSum));

    %temp jac
    ii(end+1)=ix; jj(end+1)=lastid; vv(end+1)=-ikl*(-u(ixlr)*ci(ikl,kSum)-sgn(ikl)*u(ixli)*cr(ikl,kSum));
    ii(end+1)=iy; jj(end+1)=lastid; vv(end+1)=-ikl*(-u(iylr)*ci(ikl,kSum)-sgn(ikl)*u(iyli)*cr(ikl,kSum));
    ii(end+1)=iz; jj(end+1)=lastid; vv(end+1)=-ikl*(-u(izlr)*ci(ikl,kSum)-sgn(ikl)*u(izli)*cr(ikl,kSum));
 
    %jac
    ii(end+1)=ix; jj(end+1)=ixlr; vv(end+1)=-s*(cr(ikl,kSum))               -ikl*om*(-ci(ikl,kSum));
    ii(end+1)=ix; jj(end+1)=ixli; vv(end+1)=-s*(-ci(ikl,kSum)*sgn(ikl))    -ikl*om*(-sgn(ikl)*cr(ikl,kSum));
    ii(end+1)=ix; jj(end+1)=iylr; vv(end+1)=s*(cr(ikl,kSum));
    ii(end+1)=ix; jj(end+1)=iyli; vv(end+1)=s*(-ci(ikl,kSum)*sgn(ikl));

    ii(end+1)=iy; jj(end+1)=iylr; vv(end+1)=-1*(cr(ikl,kSum))               -ikl*om*(-ci(ikl,kSum));
    ii(end+1)=iy; jj(end+1)=iyli; vv(end+1)=-1*(-ci(ikl,kSum)*sgn(ikl))    -ikl*om*(-sgn(ikl)*cr(ikl,kSum));
    ii(end+1)=iy; jj(end+1)=ixlr; vv(end+1)=r*(cr(ikl,kSum));
    ii(end+1)=iy; jj(end+1)=ixli; vv(end+1)=r*(-ci(ikl,kSum)*sgn(ikl));

    ii(end+1)=iz; jj(end+1)=izlr; vv(end+1)=-b*(cr(ikl,kSum))                 -ikl*om*(-ci(ikl,kSum));
    ii(end+1)=iz; jj(end+1)=izli; vv(end+1)=-b*(-ci(ikl,kSum)*sgn(ikl))      -ikl*om*(-sgn(ikl)*cr(ikl,kSum));


else
    g(ix)=g(ix)+s*(ci(ikl,kSum)*u(iylr)+cr(ikl,kSum)*sgn(ikl)*u(iyli));
    g(ix)=g(ix)-s*(ci(ikl,kSum)*u(ixlr)+cr(ikl,kSum)*sgn(ikl)*u(ixli));
    g(ix)=g(ix)-ikl*om*(+u(ixlr)*cr(ikl,kSum)-sgn(ikl)*u(ixli)*ci(ikl,kSum));

    g(iy)=g(iy)+r*(ci(ikl,kSum)*u(ixlr)+cr(ikl,kSum)*sgn(ikl)*u(ixli));
    g(iy)=g(iy)-1*(ci(ikl,kSum)*u(iylr)+cr(ikl,kSum)*sgn(ikl)*u(iyli));
    g(iy)=g(iy)-ikl*om*(+u(iylr)*cr(ikl,kSum)-sgn(ikl)*u(iyli)*ci(ikl,kSum));

    g(iz)=g(iz)-b*(ci(ikl,kSum)*u(izlr)+cr(ikl,kSum)*sgn(ikl)*u(izli));
    g(iz)=g(iz)-ikl*om*(+u(izlr)*cr(ikl,kSum)-sgn(ikl)*u(izli)*ci(ikl,kSum));

    %temp jac
    ii(end+1)=ix; jj(end+1)=lastid; vv(end+1)=-ikl*(+u(ixlr)*cr(ikl,kSum)-sgn(ikl)*u(ixli)*ci(ikl,kSum));
    ii(end+1)=iy; jj(end+1)=lastid; vv(end+1)=-ikl*(+u(iylr)*cr(ikl,kSum)-sgn(ikl)*u(iyli)*ci(ikl,kSum));
    ii(end+1)=iz; jj(end+1)=lastid; vv(end+1)=-ikl*(+u(izlr)*cr(ikl,kSum)-sgn(ikl)*u(izli)*ci(ikl,kSum));
 
    %jac
    ii(end+1)=ix; jj(end+1)=ixlr; vv(end+1)=-s*(ci(ikl,kSum))               -ikl*om*(cr(ikl,kSum));
    ii(end+1)=ix; jj(end+1)=ixli; vv(end+1)=-s*(cr(ikl,kSum)*sgn(ikl))     -ikl*om*(-sgn(ikl)*ci(ikl,kSum));
    ii(end+1)=ix; jj(end+1)=iylr; vv(end+1)=s*(ci(ikl,kSum));
    ii(end+1)=ix; jj(end+1)=iyli; vv(end+1)=s*(cr(ikl,kSum)*sgn(ikl));

    ii(end+1)=iy; jj(end+1)=iylr; vv(end+1)=-1*(ci(ikl,kSum))              -ikl*om*(cr(ikl,kSum));
    ii(end+1)=iy; jj(end+1)=iyli; vv(end+1)=-1*(cr(ikl,kSum)*sgn(ikl))    -ikl*om*(-sgn(ikl)*ci(ikl,kSum));
    ii(end+1)=iy; jj(end+1)=ixlr; vv(end+1)=+r*(ci(ikl,kSum));
    ii(end+1)=iy; jj(end+1)=ixli; vv(end+1)=+r*(cr(ikl,kSum)*sgn(ikl));

    ii(end+1)=iz; jj(end+1)=izlr; vv(end+1)=-b*(ci(ikl,kSum))                 -ikl*om*(cr(ikl,kSum));
    ii(end+1)=iz; jj(end+1)=izli; vv(end+1)=-b*(cr(ikl,kSum)*sgn(ikl))       -ikl*om*(-sgn(ikl)*ci(ikl,kSum));
end

end
    
end

%last eq
g(lastid)=sum(u(neq*nt+1:neq:end-1).*[0:1:nt-1]');
aa=neq*nt+1:neq:lastid-1;
ii=[ii,aa./aa*lastid]; jj=[jj,aa]; vv=[vv,[0:1:nt-1]];

jac=sparse(ii,jj,vv);
end