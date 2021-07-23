from consts import *


def phi1(deltaS, yVec):
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
    :return: vector after one exp(linear)  step
    '''
    tmp1 = phi1(1 / 2 * dt, yVec)
    tmp2 = zeta2(dt, tmp1)
    tmp3 = phi1(1 / 2 * dt, tmp2)
    return tmp3


def f(deltaTau, yVec):
    '''

    :param deltaTau: time step
    :param yVec: y vector
    :return: y mapped by f
    '''
    retVec = []
    # N must be even
    for n in range(0, N - 1, 2):
        retVec.append(
            yVec[n] * np.exp(1j * np.conj(yVec[n]) * yVec[n + 1] * deltaTau)

        )
        retVec.append(
            yVec[n + 1] * np.exp(1j * np.conj(yVec[n + 1]) * yVec[n] * deltaTau)
        )
    return retVec


def g(deltaTau, yVec):
    '''

    :param deltaTau: time step
    :param yVec:  y vector
    :return: y mapped by g
    '''
    retVec = []
    # n=0, N-1 yn does not change
    retVec.append(yVec[0])
    for n in range(1, N - 2, 2):
        retVec.append(
            yVec[n] * np.exp(1j * np.conj(yVec[n]) * yVec[n + 1] * deltaTau)
        )
        retVec.append(
            yVec[n + 1] * np.exp(1j * np.conj(yVec[n + 1]) * yVec[n] * deltaTau)
        )
    retVec.append(yVec[N - 1])

    return retVec


def NP(deltaS, yVec):
    '''

    :param deltaS: time step
    :param yVec:  y vector
    :return: vector after nonlinear step
    '''
    tmp1 = f(1 / 2 * deltaS, yVec)
    tmp2 = g(deltaS, tmp1)
    tmp3 = f(1 / 2 * deltaS, tmp2)
    return tmp3


def oneStepS2(inVec):
    '''

    :param inVec: input wvfcnt
    :return: wvfcnt after one step S2
    '''
    tmp1 = NP(1 / 2 * dt, inVec)
    tmp2 = expmidtF(tmp1)
    tmp3 = NP(1 / 2 * dt, tmp2)
    return tmp3


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