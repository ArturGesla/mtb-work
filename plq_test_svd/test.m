Q=[1,2,3;
    2,3,1;
    3,1,2];
L=diag([1,2,0])
A=(Q)*L*inv(Q)

 [U,S,V]=svd(A)
% [evc,evs]=eig(A'*A)
[evc,evs]=eig(A)
%%

Q=[1,2,3,4;
    2,3,1;
    3,1,2];
L=diag([1,2,0])
A=(Q)*L*inv(Q)

 [U,S,V]=svd(A)
% [evc,evs]=eig(A'*A)
[evc,evs]=eig(A)