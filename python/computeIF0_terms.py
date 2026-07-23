import numpy as np
from scipy.io import netcdf_file


booz_file="dat/boozmn_VMECoutput.nc"

ntheta=128
nzeta=128


f=netcdf_file(booz_file,'r')


rmnc=f.variables['rmnc_b'].data.copy()
zmns=f.variables['zmns_b'].data.copy()

xm=f.variables['ixm_b'].data.copy()
xn=f.variables['ixn_b'].data.copy()

ns=rmnc.shape[0]
s=ns-1


theta=np.linspace(0,2*np.pi,ntheta,endpoint=False)
zeta=np.linspace(0,2*np.pi,nzeta,endpoint=False)

theta2d,zeta2d=np.meshgrid(theta,zeta,indexing='ij')


# radial derivative via finite difference

ds=1

rmnc_s=(rmnc[s]-rmnc[s-1])
zmns_s=(zmns[s]-zmns[s-1])


R=np.zeros_like(theta2d)
Z=np.zeros_like(theta2d)

Rs=np.zeros_like(theta2d)
Zs=np.zeros_like(theta2d)

Rt=np.zeros_like(theta2d)
Rz=np.zeros_like(theta2d)

Zt=np.zeros_like(theta2d)
Zz=np.zeros_like(theta2d)


for m in range(len(xm)):

 angle=xm[m]*theta2d-xn[m]*zeta2d

 cos=np.cos(angle)
 sin=np.sin(angle)

 R+=rmnc[s,m]*cos
 Z+=zmns[s,m]*sin

 Rs+=rmnc_s[m]*cos
 Zs+=zmns_s[m]*sin

 Rt+=-xm[m]*rmnc[s,m]*sin
 Rz+=xn[m]*rmnc[s,m]*sin

 Zt+=xm[m]*zmns[s,m]*cos
 Zz+=-xn[m]*zmns[s,m]*cos



X=R*np.cos(zeta2d)
Y=R*np.sin(zeta2d)

Xs=Rs*np.cos(zeta2d)
Ys=Rs*np.sin(zeta2d)
Zs=Zs

Xt=Rt*np.cos(zeta2d)
Yt=Rt*np.sin(zeta2d)
Zt=Zt

Xz=Rz*np.cos(zeta2d)-R*np.sin(zeta2d)
Yz=Rz*np.sin(zeta2d)+R*np.cos(zeta2d)
Zz=Zz


es=np.stack((Xs,Ys,Zs),axis=2)
et=np.stack((Xt,Yt,Zt),axis=2)
ez=np.stack((Xz,Yz,Zz),axis=2)


g_cov=np.zeros((ntheta,nzeta,3,3))

g_cov[:,:,0,0]=np.sum(es*es,axis=2)
g_cov[:,:,0,1]=np.sum(es*et,axis=2)
g_cov[:,:,0,2]=np.sum(es*ez,axis=2)

g_cov[:,:,1,0]=g_cov[:,:,0,1]
g_cov[:,:,1,1]=np.sum(et*et,axis=2)
g_cov[:,:,1,2]=np.sum(et*ez,axis=2)

g_cov[:,:,2,0]=g_cov[:,:,0,2]
g_cov[:,:,2,1]=g_cov[:,:,1,2]
g_cov[:,:,2,2]=np.sum(ez*ez,axis=2)


g_contra=np.linalg.inv(g_cov)


h=g_contra[:,:,0,0]
tau=g_contra[:,:,0,1]
xi=g_contra[:,:,0,2]
y=g_contra[:,:,1,2]


Jac=1/np.sqrt(np.linalg.det(g_contra))


termI=(y-tau*xi)/h**2
termF=(xi**2)/h**2


weight=1/Jac


avgI=np.sum(termI*weight)/np.sum(weight)
avgF=np.sum(termF*weight)/np.sum(weight)


with open("dat/I_F0_out.dat","w") as out:

 out.write("{:.15e}\n".format(avgI))
 out.write("{:.15e}\n".format(avgF))


print("Written to I_F0_out.dat")
