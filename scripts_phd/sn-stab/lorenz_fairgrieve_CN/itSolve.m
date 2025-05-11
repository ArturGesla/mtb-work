a=bicgstab(J,g,1e-15,60);
%%
[L,U] = ilu(sparse(J),struct('type','ilutp','droptol',1e-16));
a=bicgstab(J,g,1e-15,6,L,U);
norm(a-J\g)