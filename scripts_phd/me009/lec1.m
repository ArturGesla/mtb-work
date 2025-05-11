cd /people/gesla/Documents/git/mtb-work/me009;

%%
A=[2,1,1;
    1,2,0;
    1,0,2];

[evc,evs]=eig(A);

v=evc(:,1)+evc(:,2); %rand(3,1);

for i=1:200
    u=A*v./norm(v);
    i
    norm(u)/norm(v)
    u=v;
end
%%
n=4;
A=eye(n-1)*2-diag(ones(n-2,1),-1)-diag(ones(n-2,1),1);
e1=min(eig(A))*n^2
pi-sqrt(e1)

n=8;
A=eye(n-1)*2-diag(ones(n-2,1),-1)-diag(ones(n-2,1),1);
e2=min(eig(A))*n^2
pi-sqrt(e2)

n=16;
A=eye(n-1)*2-diag(ones(n-2,1),-1)-diag(ones(n-2,1),1);
e3=min(eig(A))*n^2
% 
% % (e2-e1)/(e3-e2)
% % (sqrt(e2)-sqrt(e1))/(sqrt(e3)-sqrt(e2))
(pi-sqrt(e1))/(pi-sqrt(e2))

