import numpy as np
from math import floor
import matplotlib.pyplot as plt
from decimal import Decimal
import time

origTime = time.time()
FLOPS = 0


def firstString():
    for j in range(0, Nx - 1):
        nullMatrix[0, j] = B1
    nullMatrix[0, Nx - 1] = B2


def Borders():
    for j in range(0, Nt):
        nullMatrix[j, 0] = B1
        nullMatrix[j, Nx - 1] = B2


def nextAlpha(preAlpha):
    result = sigma / (1 + (2 * sigma) - (sigma * preAlpha))
    return result


def nextBetta(preBetta, alpha, w):
    a = 4
    result = (alpha / sigma) * (w + sigma * preBetta - (dt * (w * (1 - w) * (a - w))))
    return result


def ourFunc(alpha, betta, w):
    result = alpha * w + betta
    return result


def method(FLOPS):
    for i in range(1, Nt - 1):
        Alphas = np.zeros(Nx)
        Alphas[1] = 0
        for j in range(1, Nx - 1):
            Alphas[j + 1] = nextAlpha(Alphas[j])
            FLOPS += 5
        Bettas = np.zeros(Nx)
        Bettas[1] = B1
        for j in range(1, Nx - 1):
            Bettas[j + 1] = nextBetta(Bettas[j], Alphas[j + 1], nullMatrix[i - 1, j])
            FLOPS += 28
        for j in range(Nx - 2, 0, -1):
            nullMatrix[i + 1, j] = ourFunc(Alphas[j + 1], Bettas[j + 1], nullMatrix[i + 1, j + 1])
            FLOPS += 2
    print('FLOPS на данном шаге =', FLOPS)


B1 = 10
B2 = 0
dx = 0.5
dt = 0.01
matrixToCompare = []

pointToCompare = 11
for j in range(0, pointToCompare):
    print('dx =', dx)
    ltX = 5
    ltT = 5
    Nt = floor(ltT / dt)
    Nx = floor(ltX / dx)
    sigma = dt / (dx ** 2)
    nullMatrix = np.zeros((Nt, Nx))
    firstString()
    Borders()
    method(FLOPS)
    dx = dx/2
    matrixToCompare.append(nullMatrix[Nt - 1])
    print('Общее время расчёта =', time.time() - origTime, 'секунд')

for i in range (0, pointToCompare - 1):
    F = np.zeros(len(matrixToCompare[i]))
    for j in range (0, len(matrixToCompare[i])-1):
        F[j] = abs(matrixToCompare[i][j] - matrixToCompare[i+1][j*2])
    print('локальная ошибка при шаге', dx, '=', max(F))
endTime = time.time() - origTime
print('Общее время расчёта =', endTime, 'секунд')
print('Количество шагов, необходимое для достижения нужной точности:', pointToCompare)
print('Шаг, необходимый для достижения необходимой точности:', dx,'при dt =', dt)
J = []
for j in range(0, len(nullMatrix[0])):
    J.append(j)
plt.plot(J, nullMatrix[Nt - 1])
plt.show()
