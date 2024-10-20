function [g,jac]=calculateRhsAndJac(neq,nt,u,mu,gm,om1)
%%
% Noack system
% f=[mu*u-v*g-w*u;
%     mu*v+u/g-v*w;
%     -w+u^2+v^2*g*g]; 
%4 nonlinear terms - troche zabawy
% 

% mu=0.04; gm=1; 

ii=[]; jj=[]; vv=[];

%u=[xr0,yr0,zr0,xr1,...,xi0,yi0,zi0,...]

%Lorenz
% s=10;
% b=8/3;
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
%     nlxz=0;
%     nlxy=0;

    nlxz=0;
    nlyz=0;
    nlyy=0;
    nlxx=0;



    kSum=mod(ik,nt);
    for ikl=kSum-nt+1:nt-1
        ikr=kSum-ikl;
        
        ixrr=abs(ikr)*neq+1; ixri=(abs(ikr)+nt)*neq+1; ixlr=abs(ikl)*neq+1; ixli=(abs(ikl)+nt)*neq+1;
        iyrr=abs(ikr)*neq+2; iyri=(abs(ikr)+nt)*neq+2; iylr=abs(ikl)*neq+2; iyli=(abs(ikl)+nt)*neq+2;
        izrr=abs(ikr)*neq+3; izri=(abs(ikr)+nt)*neq+3; izlr=abs(ikl)*neq+3; izli=(abs(ikl)+nt)*neq+3;
        
        if(realEq==1)
%             nlxz=nlxz+u(ixrr)*u(izlr)- sign(ikr)* sign(ikl)*u(ixli)*u(izri);
%             nlxy=nlxy+u(ixrr)*u(iylr)- sign(ikr)* sign(ikl)*u(ixli)*u(iyri);
%             
%             %jac
%             ii(end+1)=iy; jj(end+1)=ixrr; vv(end+1)=u(izlr)*-1;
%             ii(end+1)=iy; jj(end+1)=izlr; vv(end+1)=u(ixrr)*-1;
%             ii(end+1)=iy; jj(end+1)=ixli; vv(end+1)=- sign(ikr)* sign(ikl)*u(izri)*-1;
%             ii(end+1)=iy; jj(end+1)=izri; vv(end+1)=- sign(ikr)* sign(ikl)*u(ixli)*-1;
%             
%             
%             ii(end+1)=iz; jj(end+1)=ixrr; vv(end+1)=1*u(iylr);
%             ii(end+1)=iz; jj(end+1)=iylr; vv(end+1)=u(ixrr);
%             ii(end+1)=iz; jj(end+1)=ixli; vv(end+1)=- sign(ikr)* sign(ikl)*u(iyri);
%             ii(end+1)=iz; jj(end+1)=iyri; vv(end+1)=- sign(ikr)* sign(ikl)*u(ixli);


            nlxz=nlxz+u(ixrr)*u(izlr)- sign(ikr)* sign(ikl)*u(ixli)*u(izri);
            nlyz=nlyz+u(iyrr)*u(izlr)- sign(ikr)* sign(ikl)*u(iyli)*u(izri);
            nlxx=nlxx+u(ixrr)*u(ixlr)- sign(ikr)* sign(ikl)*u(ixli)*u(ixri);
            nlyy=nlyy+u(iyrr)*u(iylr)- sign(ikr)* sign(ikl)*u(iyli)*u(iyri);

            %jac
            ii(end+1)=ix; jj(end+1)=ixrr; vv(end+1)=u(izlr)*-1;
            ii(end+1)=ix; jj(end+1)=izlr; vv(end+1)=u(ixrr)*-1;
            ii(end+1)=ix; jj(end+1)=ixli; vv(end+1)=- sign(ikr)* sign(ikl)*u(izri)*-1;
            ii(end+1)=ix; jj(end+1)=izri; vv(end+1)=- sign(ikr)* sign(ikl)*u(ixli)*-1;

            ii(end+1)=iy; jj(end+1)=iyrr; vv(end+1)=1*u(izlr)*-1;
            ii(end+1)=iy; jj(end+1)=izlr; vv(end+1)=u(iyrr)*-1;
            ii(end+1)=iy; jj(end+1)=iyli; vv(end+1)=- sign(ikr)* sign(ikl)*u(izri)*-1;
            ii(end+1)=iy; jj(end+1)=izri; vv(end+1)=- sign(ikr)* sign(ikl)*u(iyli)*-1;

            ii(end+1)=iz; jj(end+1)=ixrr; vv(end+1)=1*u(ixlr);
            ii(end+1)=iz; jj(end+1)=ixlr; vv(end+1)=u(ixrr);
            ii(end+1)=iz; jj(end+1)=ixli; vv(end+1)=- sign(ikr)* sign(ikl)*u(ixri);
            ii(end+1)=iz; jj(end+1)=ixri; vv(end+1)=- sign(ikr)* sign(ikl)*u(ixli);

            ii(end+1)=iz; jj(end+1)=iyrr; vv(end+1)=1*u(iylr)*gm*gm;
            ii(end+1)=iz; jj(end+1)=iylr; vv(end+1)=u(iyrr)*gm*gm;
            ii(end+1)=iz; jj(end+1)=iyli; vv(end+1)=- sign(ikr)* sign(ikl)*u(iyri)*gm*gm;
            ii(end+1)=iz; jj(end+1)=iyri; vv(end+1)=- sign(ikr)* sign(ikl)*u(iyli)*gm*gm;
            

            
            
        else
%             nlxz=nlxz+u(ixrr)*u(izli)* sign(ikl)+ sign(ikr)*u(ixri)*u(izlr);
%             nlxy=nlxy+u(ixrr)*u(iyli)* sign(ikl)+ sign(ikr)*u(ixri)*u(iylr);
%             
%             %jac
%             ii(end+1)=iy; jj(end+1)=ixrr; vv(end+1)=u(izli)* sign(ikl)*-1;
%             ii(end+1)=iy; jj(end+1)=izli; vv(end+1)=u(ixrr)* sign(ikl)*-1;
%             ii(end+1)=iy; jj(end+1)=ixri; vv(end+1)= sign(ikr)*u(izlr)*-1;
%             ii(end+1)=iy; jj(end+1)=izlr; vv(end+1)=sign(ikr)*u(ixri)*-1;
%             
%             ii(end+1)=iz; jj(end+1)=ixrr; vv(end+1)=u(iyli)* sign(ikl);
%             ii(end+1)=iz; jj(end+1)=iyli; vv(end+1)=u(ixrr)* sign(ikl);
%             ii(end+1)=iz; jj(end+1)=ixri; vv(end+1)= sign(ikr)*u(iylr);
%             ii(end+1)=iz; jj(end+1)=iylr; vv(end+1)= sign(ikr)*u(ixri);

            nlxz=nlxz+u(ixrr)*u(izli)* sign(ikl)+ sign(ikr)*u(ixri)*u(izlr);
            nlyz=nlyz+u(iyrr)*u(izli)* sign(ikl)+ sign(ikr)*u(iyri)*u(izlr);
            nlxx=nlxx+u(ixrr)*u(ixli)* sign(ikl)+ sign(ikr)*u(ixri)*u(ixlr);
            nlyy=nlyy+u(iyrr)*u(iyli)* sign(ikl)+ sign(ikr)*u(iyri)*u(iylr);

            %jac
            ii(end+1)=ix; jj(end+1)=ixrr; vv(end+1)=u(izli)* sign(ikl)*-1;
            ii(end+1)=ix; jj(end+1)=izli; vv(end+1)=u(ixrr)* sign(ikl)*-1;
            ii(end+1)=ix; jj(end+1)=ixri; vv(end+1)= sign(ikr)*u(izlr)*-1;
            ii(end+1)=ix; jj(end+1)=izlr; vv(end+1)=sign(ikr)*u(ixri)*-1;

            ii(end+1)=iy; jj(end+1)=iyrr; vv(end+1)=u(izli)* sign(ikl)*-1;
            ii(end+1)=iy; jj(end+1)=izli; vv(end+1)=u(iyrr)* sign(ikl)*-1;
            ii(end+1)=iy; jj(end+1)=iyri; vv(end+1)= sign(ikr)*u(izlr)*-1;
            ii(end+1)=iy; jj(end+1)=izlr; vv(end+1)= sign(ikr)*u(iyri)*-1;

            ii(end+1)=iz; jj(end+1)=ixrr; vv(end+1)=u(ixli)* sign(ikl);
            ii(end+1)=iz; jj(end+1)=ixli; vv(end+1)=u(ixrr)* sign(ikl);
            ii(end+1)=iz; jj(end+1)=ixri; vv(end+1)= sign(ikr)*u(ixlr);
            ii(end+1)=iz; jj(end+1)=ixlr; vv(end+1)= sign(ikr)*u(ixri);

            ii(end+1)=iz; jj(end+1)=iyrr; vv(end+1)=u(iyli)* sign(ikl)*gm*gm;
            ii(end+1)=iz; jj(end+1)=iyli; vv(end+1)=u(iyrr)* sign(ikl)*gm*gm;
            ii(end+1)=iz; jj(end+1)=iyri; vv(end+1)= sign(ikr)*u(iylr)*gm*gm;
            ii(end+1)=iz; jj(end+1)=iylr; vv(end+1)= sign(ikr)*u(iyri)*gm*gm;
            
        end
    end
    
    
    
    
    ikl=kSum;
    %     ixrr=abs(ikr)*neq+1; ixri=(abs(ikr)+nt)*neq+1;
    ixlr=abs(ikl)*neq+1; ixli=(abs(ikl)+nt)*neq+1;
    %         iyrr=abs(ikr)*neq+2; iyri=(abs(ikr)+nt)*neq+2;
    iylr=abs(ikl)*neq+2; iyli=(abs(ikl)+nt)*neq+2;
    %     izrr=abs(ikr)*neq+3; izri=(abs(ikr)+nt)*neq+3;
    izlr=abs(ikl)*neq+3; izli=(abs(ikl)+nt)*neq+3;
%     
%     g(ix)=s*(u(iy)-u(ix))-(realEq*(-om*kSum*u(ixli))+~realEq*(om*kSum*u(ixlr)));
%     g(iy)=r*u(ix)-u(iy)-nlxz-(realEq*(-om*kSum*u(iyli))+~realEq*(om*kSum*u(iylr)));
%     g(iz)=nlxy-b*u(iz)-(realEq*(-om*kSum*u(izli))+~realEq*(om*kSum*u(izlr)));
%     
%     %temp jac
%     ii(end+1)=ix; jj(end+1)=lastid; vv(end+1)=-(realEq*(-kSum*u(ixli))+~realEq*(kSum*u(ixlr)));
%     ii(end+1)=iy; jj(end+1)=lastid; vv(end+1)=-(realEq*(-kSum*u(iyli))+~realEq*(kSum*u(iylr)));
%     ii(end+1)=iz; jj(end+1)=lastid; vv(end+1)=-(realEq*(-kSum*u(izli))+~realEq*(kSum*u(izlr)));
%     
%     %jac
%     ii(end+1)=ix; jj(end+1)=ix; vv(end+1)=-s;
%     ii(end+1)=ix; jj(end+1)=iy; vv(end+1)=s;
%     ii(end+1)=ix; jj(end+1)=ixli; vv(end+1)=-(realEq*(-om*kSum));
%     ii(end+1)=ix; jj(end+1)=ixlr; vv(end+1)=-(~realEq*(om*kSum));
%     
%     ii(end+1)=iy; jj(end+1)=ix; vv(end+1)=r;
%     ii(end+1)=iy; jj(end+1)=iy; vv(end+1)=-1;
%     ii(end+1)=iy; jj(end+1)=iyli; vv(end+1)=-(realEq*(-om*kSum));
%     ii(end+1)=iy; jj(end+1)=iylr; vv(end+1)=-(~realEq*(om*kSum));
%     
%     ii(end+1)=iz; jj(end+1)=iz; vv(end+1)=-b;
%     ii(end+1)=iz; jj(end+1)=izli; vv(end+1)= -(realEq*(-om*kSum));
%     ii(end+1)=iz; jj(end+1)=izlr; vv(end+1)= -(~realEq*(om*kSum));


    g(ix)= mu*u(ix)-gm*u(iy)-nlxz-(realEq*(-om*kSum*u(ixli))+~realEq*(om*kSum*u(ixlr)))/om1;
    g(iy)=mu*u(iy)+u(ix)/gm-nlyz-(realEq*(-om*kSum*u(iyli))+~realEq*(om*kSum*u(iylr)))/om1;
    g(iz)=-u(iz)+nlxx+nlyy*gm*gm-(realEq*(-om*kSum*u(izli))+~realEq*(om*kSum*u(izlr)));


%     %temp jac
    ii(end+1)=ix; jj(end+1)=lastid; vv(end+1)=-(realEq*(-kSum*u(ixli))+~realEq*(kSum*u(ixlr)))/om1;
    ii(end+1)=iy; jj(end+1)=lastid; vv(end+1)=-(realEq*(-kSum*u(iyli))+~realEq*(kSum*u(iylr)))/om1;
    ii(end+1)=iz; jj(end+1)=lastid; vv(end+1)=-(realEq*(-kSum*u(izli))+~realEq*(kSum*u(izlr)));
%     
%     %jac
    ii(end+1)=ix; jj(end+1)=ix; vv(end+1)=mu;
    ii(end+1)=ix; jj(end+1)=iy; vv(end+1)=-gm;
    ii(end+1)=ix; jj(end+1)=ixli; vv(end+1)=-(realEq*(-om*kSum))/om1;
    ii(end+1)=ix; jj(end+1)=ixlr; vv(end+1)=-(~realEq*(om*kSum))/om1;

    ii(end+1)=iy; jj(end+1)=ix; vv(end+1)=1/gm;
    ii(end+1)=iy; jj(end+1)=iy; vv(end+1)=mu;
    ii(end+1)=iy; jj(end+1)=iyli; vv(end+1)=-(realEq*(-om*kSum))/om1;
    ii(end+1)=iy; jj(end+1)=iylr; vv(end+1)=-(~realEq*(om*kSum))/om1;

    ii(end+1)=iz; jj(end+1)=iz; vv(end+1)=-1;
    ii(end+1)=iz; jj(end+1)=izli; vv(end+1)= -(realEq*(-om*kSum));
    ii(end+1)=iz; jj(end+1)=izlr; vv(end+1)= -(~realEq*(om*kSum));


    
    
    
end

%last eq
g(lastid)=sum(u(neq*nt+1:neq:end-1).*[0:1:nt-1]');
aa=neq*nt+1:neq:lastid-1;
ii=[ii,aa./aa*lastid]; jj=[jj,aa]; vv=[vv,[0:1:nt-1]];

jac=sparse(ii,jj,vv);
end