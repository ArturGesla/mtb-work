B=A;
for i=2:length(xc)-1
    A(i,i)=-2/dx/dx;
    A(i,i+1)=1/dx/dx;
    A(i,i-1)=1/dx/dx;
%     b(i)=f(xc(i));
B(i,i)=1;
B(i,i+1)=1/10;
B(i,i-1)=1/10;
end
A=A*6/5;
fb=f(xc');

% A(1,1:2)=[1 1 ]; A(end,end-1:end)=[1 1 ];
A(1,1)=[1]; A(end,end)=[1 ];

%what a magic