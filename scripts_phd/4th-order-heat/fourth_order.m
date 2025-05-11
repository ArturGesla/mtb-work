
for i=3:length(xc)-2
    A(i,i-2)=-1/dx/dx/12;
    A(i,i-1)=16/dx/dx/12;
    A(i,i)=-30/dx/dx/12;
    A(i,i+1)=16/dx/dx/12;
    A(i,i+2)=-1/dx/dx/12;
    b(i)=f(xc(i));
end

for i=[2,length(xc)-1]
    A(i,i)=-2/dx/dx;
    A(i,i+1)=1/dx/dx;
    A(i,i-1)=1/dx/dx;
    b(i)=f(xc(i));
end

% A(1,1:2)=[1 1 ]; A(end,end-1:end)=[1 1 ];
A(1,1)=[1]; A(end,end)=[1 ];
