%%
A=sparse(J); b=g;
%%
tic; b\A; toc;
%%
close all;
[L,U]=ilu(A,struct('type','ilutp','droptol',1e-2));
tic; [x,flag,relres,iter,resvec]=gmres(A,b,[],1e-12,100,L,U); toc;
semilogy(resvec); hold on;
%%
tic; [x,flag,relres,iter,resvec]=gmres(A,b,[],0,100); toc;
semilogy(resvec); hold on;