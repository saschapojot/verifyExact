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

    retVec = []
    for n in range(0, N):
        retVec.append(
            yVec[n] * np.exp(-1j * deltaS * a * n)
        )

    return retVec

def phi21(deltaTau, yVec):
    '''

    :param deltaTau: time step
    :param yVec: y vector
    :return: y mapped by phi21
    '''

    retVec = []
    for n in range(0, N - 1):
        retVec.append(yVec[n] + 1j * deltaTau * yVec[n + 1])
    retVec.append(yVec[N - 1])
    return retVec


def phi22(deltaTau, yVec):
    '''

    :param deltaTau: time step
    :param yVec: y vector
    :return: y mapped by phi22
    '''
    retVec = []
    retVec.append(yVec[0])
    for n in range(1, N):
        retVec.append(
            yVec[n] + 1j * deltaTau * yVec[n - 1]
        )

    return retVec



def zeta2(deltaS, yVec):
    '''

    :param deltaS: time step
    :param yVec: y vector
    :return: y mapped by zeta2
    '''
    tmp1 = phi21(1 / 2 * deltaS, yVec)
    tmp2 = phi22(deltaS, tmp1)
    tmp3 = phi21(1 / 2 * deltaS, tmp2)
    return tmp3




def expmidtF(yVec):
    '''
    linear step in S2
    :param yVec:
    :return:
    '''
    tmp1 = phi1(1 / 2 * dt, yVec)
    tmp2 = zeta2(dt, tmp1)
    tmp3 = phi1(1 / 2 * dt, tmp2)
    return tmp3


def oneStepS2(inVec):
    '''

    :param inVec: input wvfcnt
    :return: wvfcnt after one step S2
    '''
    tmp1=oneStepRk2(inVec)
    tmp2=expmidtF(tmp1)
    psiqNext=oneStepRk2(tmp2)

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