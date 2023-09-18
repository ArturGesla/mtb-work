%% Fairgrieve part
% A=J;
% B(np*neq,np*neq)=1; B(np*neq-1,np*neq-1)=1;
% B(np*neq-2,np*neq-2)=1; B(np*neq-3,np*neq-3)=1;
% B(np*neq+1,np*neq+1)=0; B=sparse(B);
% [evc,evs]=eigs(A,B,neq,"smallestabs"); evs=diag(evs)
% lam=1./(1-evs) % one of lambda should be 1
% exp=1./T*log(abs(lam))
%% Fairgrieve part
A=J(1:end-1,1:end-1); B=[];
B(np*neq,np*neq)=1; B(np*neq-1,np*neq-1)=1;
B(2,np*neq-2)=1; B(1,np*neq-3)=1;
B(np*neq+1,np*neq+1)=0; B=B(1:end-1,1:end-1); %B=sparse(B);
% [evc,evs]=eigs(A,B,neq,"smallestabs"); evs=diag(evs)
% lam=1./(1-evs) % one of lambda should be 1
% exp=1./T*log(abs(lam))

%
[evc,evs]=eig(A,B); evs=diag(evs)
lam=1./(1-evs) % one of lambda should be 1
exp=1./T*log(abs(lam))