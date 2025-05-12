import numpy as np

a=np.loadtxt("c")

print(a)

vel=np.array([a[:,0]*0, a[:,0]*0,a[:,0]*0]).transpose()

print(vel)

print(len(a))

for i in range(0,len(a)):
	x=a[i,0];
	y=a[i,1];
	r=np.sqrt(x*x+y*y)
	#uth=(r*(2-r*r)*(r<1) + 1/r*(r>1))*1.095
	#uz=1
	#susan 2006
	om0=0.31765; om1=-0.62888; om2=2.2545; u0=0.30697; u1=0.01056; u2=-0.31889; r1=0.46643; r2=0.13051;

	uth=0;
	uth+=om0*r;
	uth+=om1*r1*r1/r*(1-np.exp(-r*r/r1/r1));
	uth+=om2*r2*r2/r*(1-np.exp(-r*r/r2/r2));

	uz=0;
	uz+=u0;
	uz+=u1*np.exp(-r*r/r1/r1);
	uz+=u2*np.exp(-r*r/r2/r2);

	th=np.arctan2(y,x)
	ux=-np.sin(th)*uth
	uy=np.cos(th)*uth
	vel[i,0]=ux
	vel[i,1]=uy
	vel[i,2]=uz

print(vel)
np.savetxt("vel",vel)
