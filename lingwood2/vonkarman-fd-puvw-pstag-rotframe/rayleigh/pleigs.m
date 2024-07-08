a=polyeig(jac0,jac1,jac2)
n=length(x)
%%
A=[zeros(n,n), eye(n);jac0,jac1];
B=[-eye(n),zeros(n,n);zeros(n,n),jac2];
b=eig(full(A),-full(B));
c=eigs(full(A),-full(B),20,'smallestabs');

%%
plot(b,'o'); hold on; plot(a,'x'); plot(c,'sq'); 

%%

a=polyeig(jac0,jac1,jac2,jac3);
n=length(x);
A=[zeros(n,n), eye(n),zeros(n,n);
    zeros(n,n),zeros(n,n),eye(n);
    jac0,jac1,jac2];
B=[-eye(n),zeros(n,n),zeros(n,n);
    zeros(n,n),-eye(n),zeros(n,n);
    zeros(n,n),zeros(n,n),jac3];


b=eig(full(A),-full(B));
c=eigs(full(A),-full(B),100,'smallestabs');

%
% a=polyeig(jac0,jac1,jac2,jac3);
% n=length(x);
% A=[zeros(n,n), eye(n),zeros(n,n);
%     zeros(n,n),zeros(n,n),eye(n);
%     jac0,jac1,zeros(n,n)];
% B=[-eye(n),zeros(n,n),zeros(n,n);
%     zeros(n,n),-eye(n),zeros(n,n);
%     zeros(n,n),zeros(n,n),jac2];
% 
% 
% b=eig(full(A),-full(B));
% c=eigs(full(A),-full(B),20,'smallestabs');

%%
clf;

a=polyeig(jac0,jac1,jac2,jac3);
[b1,b2]=polyeigs(jac0,jac1,jac2,jac3,100,0.1); b=diag(b2);
plot(b,'o'); hold on; plot(a,'x'); 
% plot(c,'sq'); 

% xlim([0 1]); ylim([-1 1])

