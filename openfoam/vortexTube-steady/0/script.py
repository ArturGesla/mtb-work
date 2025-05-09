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
	uth=r*(2-r*r)*(r<1) + 1/r*(r>1)
	th=np.arctan2(x,y)
	ux=-np.sin(th)*uth
	uy=np.cos(th)*uth
	vel[i,0]=ux
	vel[i,1]=uy
	vel[i,2]=1

print(vel)
np.savetxt("vel",vel)
