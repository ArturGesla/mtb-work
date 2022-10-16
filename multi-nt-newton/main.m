neq=2; mu=-0.015; b=-10;omega=1;

g=zeros(2,1);
J=zeros(2);
u=zeros(2,1)+0.1;

%% calc J and g
for i=1:40
    x=u(1); y=u(2); r=sqrt(x^2+y^2);
    g(1)=(mu+r^2-r^4)*x-y*(omega+b*r^2);
    g(2)=(mu+r^2-r^4)*y+x*(omega+b*r^2);
    
    J(1,1)=(mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x;
    J(1,2)=x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y);
    
    J(2,1)=y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x);
    J(2,2)=(mu+r^2-r^4)+y*(2*y-2*r^2*2*y)+x*b*2*y;
    
    du=-J\g;
    norm(du)
    u=u+du;
end