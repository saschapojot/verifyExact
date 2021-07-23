from consts import *


def zerothVector(inVec):
    '''

    :param inVec: input wvfcnt
    :return: K0 vector
    '''
    K0 = []
    K0.append(
        1j * np.abs(inVec[0]) ** 2 * inVec[1]
    )
    # n=1,...,N-2
    for n in range(1, N - 1):
        K0.append(
            1j * np.abs(inVec[n]) ** 2 * (inVec[n - 1] + inVec[n + 1])
        )
    K0.append(
        1j * np.abs(inVec[N - 1]) ** 2 * inVec[N - 2]
    )
    return K0


def firstVector(inVec, K0):
    '''

    :param inVec: input wvfcnt
    :param K0: K0 vector
    :return: K1 vector
    '''

    K1 = []

    K1.append(
        1j * np.abs(inVec[0] + 1 / 2 * dt * K0[0]) ** 2
        * (inVec[1] + 1 / 2 * dt * K0[1])
    )
    # n=1,...,N-2
    for n in range(1, N - 1):
        K1.append(
            1j * np.abs(inVec[n] + 1 / 2 * dt * K0[n]) ** 2
            * (inVec[n - 1] + 1 / 2 * dt * K0[n - 1] + inVec[n + 1] + 1 / 2 * dt * K0[n + 1])
        )

    K1.append(
        1j * np.abs(inVec[N - 1] + 1 / 2 * dt * K0[N - 1]) ** 2
        * (inVec[N - 2] + 1 / 2 * dt * K0[N - 2])
    )


    return K1



def oneStepRk2(inVec):
    '''

    :param inVec: input wvfcnt
    :return: wvfcnt after one step rk2
    '''
    K0=zerothVector(inVec)
    K1=firstVector(inVec,K0)

    out=[]
    for n in range(0,N):
        out.append(
            inVec[n]+1/4*dt*(K0[n]+K1[n])
        )

    return out



def phi1(deltaS,yVec):
    '''

    :param deltaS: time step
    :param yVec: y vector
    :return: y mapped by phi1
    '''

    for n in range(0,N):
        yVec[n]*=np.exp(-1j*deltaS*a*n)

def phi21(deltaTau, yVec):
    '''

    :param deltaTau: time step
    :param yVec: y vector
    :return: y mapped by phi21
    '''
    for n in range(0,N-1):
        yVec[n]+=1j*deltaTau*yVec[n+1]


def phi22(deltaTau, yVec):
    '''

    :param deltaTau: time step
    :param yVec: y vector
    :return: y mapped by phi22
    '''
    for n in range(1,N):
        yVec[n]+=1j*deltaTau*yVec[n-1]



def zeta2(deltaS, yVec):
    '''

    :param deltaS: time step
    :param yVec: y vector
    :return: y mapped by zeta2
    '''
    phi21(1/2*deltaS,yVec)
    phi22(deltaS,yVec)
    phi21(1/2*deltaS,yVec)




def expmidtF(yVec):
    '''
    linear step in S2
    :param yVec:
    :return:
    '''
    phi1(1/2*dt,yVec)
    zeta2(dt,yVec)
    phi1(1/2*dt,yVec)


def oneStepS2(inVec):
    '''

    :param inVec: input wvfcnt
    :return: wvfcnt after one step S2
    '''
    psi1q=oneStepRk2(inVec)
    expmidtF(psi1q)
    psiqNext=oneStepRk2(psi1q)

    return psiqNext


def l2Dist(aVec,bVec):
    '''

    :param aVec: wvfnction 1
    :param bVec: wvfunction 2
    :return: l2 distance of aVec and bVec
    '''

    n2=0
    for n in range(0,N):
        n2+=np.abs(aVec[n]-bVec[n])**2

    return np.sqrt(n2)


def lInfDist(aVec,bVec):
    '''

    :param aVec: wvfnction 1
    :param bVec: wvfunction 2
    :return: linf distance of aVec and bVec
    '''

    diffAbs=[np.abs(aVec[n]-bVec[n]) for n in range(0,N)]
    return np.max(diffAbs)