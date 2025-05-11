function [f] = reshapeplot(X,geo,prob,nbfig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    f.U = reshape(X(1:prob.nvar:end-3),geo.nx,geo.ny);
    f.V = reshape(X(2:prob.nvar:end-2),geo.nx,geo.ny);
    f.P = reshape(X(3:prob.nvar:end-1),geo.nx,geo.ny);
    f.W = reshape(X(4:prob.nvar:end),geo.nx,geo.ny);
    f.U(end,:) = [];
    f.V(:,end) = [];
    
    
    if nbfig >0
        
        figure(nbfig)
        subplot(2,2,1)
        contour(geo.xu,geo.yc,f.U.')
        title([' Ray = : ',num2str(prob.Ray),' U max : ',num2str(max(f.U(:))),' min : ',num2str(min(f.U(:))) ]);
        axis image
        subplot(2,2,2)
        contour(geo.xc,geo.yv,f.V')
        title(['V max : ',num2str(max(f.V(:))),' min : ',num2str(min(f.V(:))) ]);
        axis image
        subplot(2,2,3)
        contour(geo.xc(2:end-1),geo.yc(2:end-1),f.P(2:end-1,2:end-1).')
        title(['P max :',num2str(max(max(f.P( 2:end-1,2:end-1 )))),'min : ',num2str(min(min((f.P(2:end-1,2:end-1)) )))   ]);
        axis image
        drawnow
        subplot(2,2,4)
        contour(geo.xc(2:end-1),geo.yc(2:end-1),f.W(2:end-1,2:end-1).')
        title(['W max : ',num2str(max(max(f.W(2:end-1,2:end-1)))),'min : ',num2str(min(min((f.W(2:end-1,2:end-1)) )))   ]);
        axis image
        drawnow
        %pause
    end

end

