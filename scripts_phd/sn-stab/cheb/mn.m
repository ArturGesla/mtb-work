clc; clear; 
%%

x=0:pi*1e-4:2*pi;
y=sin(2*x); x=x/pi-1;

yn=0;
for k=0:8
    yn=yn+chebCoeff(x,y,k)*cos(k*acos(x));
end

%
plot(x,[y;yn]')

%%



u=[0 0 0 1 0 0 0 ];
n=length(u)-1;
g=u*0;
for ik=0:n
    
    ip=ik+1;
    %temporal
    for ikl=ik+1:2:n
    ipl=ikl+1;    
    a=u(ipl);
%     g(ip)=g(ip)-a*ikl*(2-(1*(ik==0)));
    end
    
    %linear
    ikl=ik;
    ipl=ikl+1;    
    
%     g(ip)=g(ip)+u(ipl);
    
    %non linear
    for ikl=ik:-1:0
        ikr=ik-ikl;
        
        ipl=ikl+1;
        ipr=ikr+1;
        
        g(ip)=g(ip)+u(ipl)*u(ipr);
        
    end
    for ikl=ik:n
        ikr=ikl-ik;
        
        ipl=ikl+1;
        ipr=ikr+1;
        
        g(ip)=g(ip)+u(ipl)*u(ipr)*(1+(ik~=0));
        
    end
    
end

%%




x = cos(pi*(0:0.02:2)/2)'; % establish 3 Chebyshev grid points
% V = exp(x); % evaluate f(x) at Chebyshev grid points
V = cos(x*pi)+sin(x*pi); % evaluate f(x) at Chebyshev grid points
%
a = fcht(V);

%

v2=[V;flipud(V(2:end-1))]; z=real(fft(v2)./length(v2)); z(2:end)=2*z(2:end);
% 

xx = linspace(-1,1); % create dense grid over domain
% g = a(1)*1 + a(2)*xx + a(3)*(2*xx.^2 - 1); % sum the first three Chebyshev
g=0;
for i=0:5
    g=g+z(i+1)*cos(i*acos(xx));
end
% polynomials with respect to their corresponding weights
plot(x,V,xx,g); % visualize the approximation