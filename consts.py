import numpy as np
import matplotlib.pyplot as plt

bt=0.4

N=2**11

a=16/N

x0=N/2
tTot=20
Q=2**14
dt=tTot/Q

def psiN(n, t):
    '''

    :param n: lattice site
    :param t: time
    :return: exact value of wvfunction
    '''

    xVal = 2 / (a * bt) * np.sinh(bt) * (np.cos(a * t) - 1) + x0
    phiVal = -2 / a * np.cosh(bt) * np.sin(a * t)
    return np.sinh(bt) * 1 / np.cosh(bt * (n - xVal)) * np.exp(-1j * (phiVal + a * n * t))

#init time
t0=0
psi0=[psiN(n,t0) for n in range(0,N)]

#
# plt.figure()
# plt.plot(range(0,N), np.abs(psi0)**2)
#
# plt.show()
# plt.close()