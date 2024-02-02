ip=30; %two neg 20 30 one one
x=u(ip*neq-2);
y=u(ip*neq-1);
z=u(ip*neq);
q=[x,y,z];

n=reshape(reshape(gstab,[neq,np])./vecnorm(reshape(gstab,[neq,np])),[neq*np,1]);
P=eye(length(gstab))-n*n';
mask=[];
for i=1:np
mask=blkdiag(mask,ones(neq,neq));
end
P=P.*mask;
A=P*Jstab*P;
A2=A(ip*neq-2:ip*neq,ip*neq-2:ip*neq);
[v,ev]=eig(A2);
close all;
utang=v(:,1);plot3([x,x+utang(1)],[y,y+utang(2)],[z,z+utang(3)],'m'); hold on;
utang=v(:,2);plot3([x,x+utang(1)],[y,y+utang(2)],[z,z+utang(3)],'c');
utang=v(:,3);plot3([x,x+utang(1)],[y,y+utang(2)],[z,z+utang(3)],'b');

%
%
% close all;
rng2=2;
plot3(u(ip*neq-2-rng2*neq:neq:ip*neq-2+rng2*neq),u(ip*neq-1-rng2*neq:neq:ip*neq-1+rng2*neq),u(ip*neq-rng2*neq:neq:ip*neq+rng2*neq));
grid on; hold on; axis equal;
% plot3([0,x],[0,y],[0,z])
utang=gstab(ip*neq-2:ip*neq);
utang=utang*0.01;
plot3([x,x+utang(1)],[y,y+utang(2)],[z,z+utang(3)],'r');
legend(num2str(ev(1,1)),num2str(ev(2,2)),num2str(ev(3,3)),"orbit","orbit tangent")
% plot3([0,utang(1)],[0,utang(2)],[0,utang(3)],'k');
view(-utang)
% view([0,0,1])
xlabel("x"); ylabel("y")

%
dx=0.1;
dy=0.1;
rng=5;
for i=-rng:rng
    for j=-rng:rng
     g=evalG(x+i*dx,y+j*dy,z,sigma,r,b);
     utang=g*0.02;
     plot3([x+i*dx,x+i*dx+utang(1)],[y+j*dy,y+j*dy+utang(2)],[z,z+utang(3)],'k-','HandleVisibility','off'); hold on;
     plot3([x+i*dx+utang(1)],[y+j*dy+utang(2)],[z+utang(3)],'b>','HandleVisibility','off'); hold on;
    end
end
        g=evalG(x,y,z,sigma,r,b);
        view(-g)
        title("Lorenz PO | r: "+num2str(r)+" sigma: "+num2str(sigma)+" b: "+num2str(b))
         
        set(gcf,"Position",[   1          41        1366         651]);
        %%
        exportgraphics(gcf,"phase5.png","Resolution",150)