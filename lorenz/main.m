r=28;
lorenz = @(t,x) [10*(x(2)-x(1));r*x(1)-1*x(1)*x(3)-x(2);x(1)*x(2)-(8/3)*x(3)];    % Anonymous Function
[T,X] = ode45(lorenz, [0 1.5586], [-6.2262, -11.0027,13.0515]);
%%
hold on;
grid on;
 plot3(X(:,1),X(:,2),X(:,3))
 plot3(X(1,1),X(1,2),X(1,3),'ro')
 disp(sqrt(8/3*(r-1)))
 %%
plot(T,X)

%% fit
x=X(:,1);
y=X(:,2);
z=X(:,3);

dx=gradient(x,T);
dy=gradient(y,T);
dz=gradient(z,T);

plot3(x,y,dx)

n=length(x);

A=[x,y,z,x.^2,x.*y,y.^2,x.*z,y.*z,z.^2];
b=dx;
cx=inv(A'*A)*(A'*b);
b=dy;
cy=inv(A'*A)*(A'*b);
b=dz;
cz=inv(A'*A)*(A'*b);

%%
disp(["1","x","y","z","xx","xy","yy","xz","yz","zz"])
disp(cx')
disp(cy')
disp(cz')

c=[cx';cy';cz'];
round(c)

