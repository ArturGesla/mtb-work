ARotorStatorUVPmodereducedntQuickHack;
clc; clear; 
%%
a=load("jacL_23.mat"); jacL=a.SA; %jacL=a.SA(1:end-1,1:end-1);
%rotst -mode 8 -nx 2 -nz 3 -restart 3000000 -contipath 1 3000000 -nModes 3 -noSolve 2

%%

a=importdata("jac.dat"); a=a.data; jacA=sparse(a(2:end,1),a(2:end,2),a(2:end,3)); jacA=jacA(1:end-1,1:end-1);
%%
 %varible order
% a1=[3;2;4;1]; a=[];
a1=[3;2;4;1]; a=[];
% for i=1:5*4*6
for i=1:4*3*2
    a=[a,a1+((i-1)*4)];
end
% a2=reshape(a,[5*4*4*6,1]);
a2=reshape(a,[4*3*2*4,1]);

%% ik order 
a=[];
a1=[0;3;1;4;2;5];
% for i=1:5*4
for i=1:4*3
    for ik=1:length(a1)
        a=[a;[1:4]'+(i-1)*4*6+a1(ik)*4];
    end
end
%%
spy(jacL(a2,a2),'b+'); hold on; spy(jacA,'bx'); spy(jacL(a2,a2)-jacA,'rsq');
%%
l2=jacL(a2,a2);

plot(l2(43,:)); hold on; plot(jacA(43,:))
% plot(l2(43,:)-jacA(43,:))


