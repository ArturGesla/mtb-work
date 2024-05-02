d=0.024
h=0.002

%to meter
d=d*1e3;
h=h*1e3;




fid=fopen('cylinder.stl','w');
fprintf(fid,'solid cylinder \n');
nt=100;
r=d/2;
phi=0:2*pi/nt:2*pi;
x=r*cos(phi); y=r*sin(phi); 

%top
for i=1: length(phi)-1
%     for i=1
z=h/2;
v0=[0;0;z];
v1=[x(i),y(i),z];
v2=[x(i+1),y(i+1),z];
a=v1-v0;
b=v2-v0;
normal=cross(a,b);

fprintf(fid,'facet normal %2.6e %2.6e %2.6e \n', normal(1),normal(2),normal(3));

fprintf(fid,'\touter loop\n');
fprintf(fid,'\t\tvertex %2.6e %2.6e %2.6e \n',v0(1),v0(2),v0(3));
fprintf(fid,'\t\tvertex %2.6e %2.6e %2.6e \n',v1(1),v1(2),v1(3));
fprintf(fid,'\t\tvertex %2.6e %2.6e %2.6e \n',v2(1),v2(2),v2(3));

fprintf(fid,'\tendloop\n');
fprintf(fid,'endfacet\n');

%side
%t1
v0=[x(i), y(i),h/2];
v1=[x(i), y(i),-h/2];
v2=[x(i+1), y(i+1),-h/2];
a=v1-v0;
b=v2-v0;
normal=cross(a,b)

fprintf(fid,'facet normal %2.6e %2.6e %2.6e \n', normal(1),normal(2),normal(3));

fprintf(fid,'\touter loop\n');
fprintf(fid,'\t\tvertex %2.6e %2.6e %2.6e \n',v0(1),v0(2),v0(3));
fprintf(fid,'\t\tvertex %2.6e %2.6e %2.6e \n',v1(1),v1(2),v1(3));
fprintf(fid,'\t\tvertex %2.6e %2.6e %2.6e \n',v2(1),v2(2),v2(3));

fprintf(fid,'\tendloop\n');
fprintf(fid,'endfacet\n');


%t2
v0=[x(i), y(i),h/2]
v1=[x(i+1), y(i+1),-h/2]
v2=[x(i+1), y(i+1),+h/2]
a=v1-v0
b=v2-v0
normal=cross(a,b)

fprintf(fid,'facet normal %2.6e %2.6e %2.6e \n', normal(1),normal(2),normal(3));

fprintf(fid,'\touter loop\n');
fprintf(fid,'\t\tvertex %2.6e %2.6e %2.6e \n',v0(1),v0(2),v0(3));
fprintf(fid,'\t\tvertex %2.6e %2.6e %2.6e \n',v1(1),v1(2),v1(3));
fprintf(fid,'\t\tvertex %2.6e %2.6e %2.6e \n',v2(1),v2(2),v2(3));

fprintf(fid,'\tendloop\n');
fprintf(fid,'endfacet\n');


%bottom
z=-h/2
v0=[0;0;z]
v2=[x(i), y(i),z]
v1=[x(i+1), y(i+1),z]
a=v1-v0
b=v2-v0
normal=cross(a,b)

fprintf(fid,'facet normal %2.6e %2.6e %2.6e \n', normal(1),normal(2),normal(3));

fprintf(fid,'\touter loop\n');
fprintf(fid,'\t\tvertex %2.6e %2.6e %2.6e \n',v0(1),v0(2),v0(3));
fprintf(fid,'\t\tvertex %2.6e %2.6e %2.6e \n',v1(1),v1(2),v1(3));
fprintf(fid,'\t\tvertex %2.6e %2.6e %2.6e \n',v2(1),v2(2),v2(3));

fprintf(fid,'\tendloop\n');
fprintf(fid,'endfacet\n');


end

fprintf(fid,'endsolid cylinder \n')
fclose(fid)
