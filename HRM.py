import numpy as np
from functools import partial
import math
from scipy import stats

sphiy = 1000.0
mphiy = 0.0
nsy2 = 0.0001
nssy2 = 0.0001


def getA(i, j, m = 5):
    if i* m <= j and j <= (i + 1) * m - 1:
        return 1 / m
    else:
        return 0

def arVar(i, j, phiY = 0.3, sigmaYSquared = 1):
    return sigmaYSquared * math.pow(phiY, i - j) / (1 - math.pow(phiY, 2))

class HRM:
    def __init__(self, m, nx, ny, gibbsIter):
        self.m = m
        self.nx = nx
        self.ny = ny
        self.A = np.fromfunction(np.vectorize(partial(getA, m = self.m)), (ny, nx))

        self.gibbsIter = gibbsIter

        self.phiYVals = [0.0]
        self.sigmaY2Vals = [0.1]
        self.lambdaVals = []



    def coarseOneStep(self, oldPhiY, oldSigmay2, yVals):
        sumy2 = 0
        sumyby = 0
        for i in yVals:
            sumy2 += i * i
        for i in range(1, len(yVals)):
            sumyby += yVals[i - 1] * yVals[i]

        phiVar = 1 / ((sumy2 - yVals[0] * yVals[0]) / oldSigmay2 + 1 / sphiy)
        phiMean = phiVar * (sumyby / oldSigmay2 + mphiy / sphiy)

        norm = stats.norm(loc=phiMean, scale=math.sqrt(phiVar))
        nextPhi = norm.rvs(1)
        while nextPhi >= 1 or nextPhi <= -1:
            nextPhi = norm.rvs(1)

        self.phiYVals.append(nextPhi)

        naux = nsy2 + len(yVals) - 1
        nsaux = nssy2 + (sumy2 - yVals[0] * yVals[0]) - 2 * oldPhiY * sumyby + oldPhiY * oldPhiY * (
        sumy2 - yVals[-1] * yVals[-1])

        sigmay2 = nsaux / np.random.chisquare(naux)
        self.sigmaY2Vals.append(sigmay2)

        return nextPhi, sigmay2

    def lambdaAndExpCoarse(self, oldPhiY, oldSigmay2):
        pass





