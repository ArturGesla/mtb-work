i=1;
%%
u=rand(np*3+1,1);
for i=1:length(u)
uorg=u;



eps=1e-6;

% i=2;
up=u; um=u;
up(i)=u(i)+eps;
um(i)=u(i)-eps;

u=up;
evalJacRhs;
gp=g;

 u=um;
evalJacRhs;
 gm=g;

evalJacRhs;
J=full(J);
% u=uorg;

col=(gp-gm)/2/eps;
col=col';

J(:,i);
i
norm(col-J(:,i))
if(norm(col-J(:,i))>1e-6) 
    break;
end

dcol=col'-J(:,i);
end
% rank(full(J))
% i=i+1