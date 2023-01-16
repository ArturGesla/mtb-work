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
