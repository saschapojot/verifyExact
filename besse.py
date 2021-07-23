import numpy as np
import matplotlib.pyplot as plt
tTot=100

Q=2**12
dt=tTot/Q
psi0=1

#init
V1=np.abs(psi0)**2
V2=2*np.abs(psi0)**2-V1
psiAll=[]
psiAll.append(psi0)
for q in range(0,Q):
    V2=2*np.abs(psiAll[q])**2-V1
    psiNext=(1+1j/2*dt   *np.abs(V2)**2)/(1-1j/2*dt*np.abs(V2)**2)*psiAll[q]
    psiAll.append(psiNext)
    V1=V2


def exact(t):
    return psi0*np.exp(1j*np.abs(psi0)**2*t)

eS=[exact(dt*q) for q in range(0,Q+1)]

l2=[np.abs(psiAll[q]-eS[q])**2 for q in range(0,Q+1)]
tAll=[dt*q for q in range(0,Q+1)]


plt.figure()
plt.plot(tAll,l2,color="black")
plt.savefig("dist.png")
plt.close()