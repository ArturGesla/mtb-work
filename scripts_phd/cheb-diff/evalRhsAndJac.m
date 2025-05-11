% function [g,jac]=evalRhsAndJac(rp,zp,ru,zu,rv,zv,rw,zw,N,R,H,u)
function evalRhsAndJac()

load input.mat;

g=u*0;
jac=0;

% %conti du/dr *2/R + u *2/(R(r+1)) + dw/dz *2/H =0
% for ix=1:length(rp)
%     for iy=1:length(zp) %grid on conti
%         xc=rp(ix); zc=zp(iy);
%         ip=((iy)+(ix)*(N+1))*4+1;
% 
%         for ikx=0:length(ru)-1 %expansion of velocity
%             for iky=0:length(zu)-1
%                 ikp=(iky+(ikx)*(N+1))*4+2;
% 
%                 g(ip)=g(ip)+u(ikp)*dTn(xc,ikx)*Tn(zc,iky)*2/R;
%                 g(ip)=g(ip)+u(ikp)*Tn(xc,ikx)*Tn(zc,iky)*2/(R*(xc+1));
%             end
%         end
% 
%         for ikx=0:length(rw)-1 %expansion of velocity
%             for iky=0:length(zw)-1
% 
%                 ikp=(iky+(ikx)*(N+1))*4+4;
% 
%                 g(ip)=g(ip)+u(ikp)*Tn(xc,ikx)*dTn(zc,iky)*2/H;
% 
%             end
%         end
% 
%         %
% 
%     end
% end

% tic;
%r mom
for ix=2:length(ru)-1
    for iy=2:length(zu)-1 %grid on r mom
        xc=ru(ix); zc=zu(iy);
        ip=((iy)+(ix)*(N+1))*4+2;

        %nl term u * du/dr *2/R + w * du/dz -v^2/r

        for ikxl=0:length(ru)-1
            for ikyl=0:length(zu)-1
                ikpl=(ikyl+(ikxl)*(N+1))*4+2;

                for ikxr=0:length(ru)-1
                    for ikyr=0:length(zu)-1
                        ikpr=(ikyr+(ikxr)*(N+1))*4+2;

                        g(ip)=g(ip)+u(ikpl)*Tn(xc,ikxl)*Tn(zc,ikyl)*u(ikpr)*dTn(xc,ikxr)*Tn(zc,ikyr)*2/R;

                    end
                end

            end
        end

    end
end
% toc;

save output.mat;



end