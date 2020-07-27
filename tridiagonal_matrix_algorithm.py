import numpy as np
from math import floor
import matplotlib.pyplot as plt
from decimal import Decimal
import time

origTime = time.time()
# Задаём ГУ для рассчёта
B1 = np.double(1)
B2 = np.double(0)
# шаги по пространству и по времени
dx = 0.05
dt = 0.0005
# длина пути, который будет разбит на шаги
ltX = 4
ltT = 4
# количество временных и пространственных узлов
Nt = floor(ltT / dt)
Nx = floor(ltX / dx)
# число Куранта, рассчитаное для сходимости схемы (1/2 >= sigma > 0)
sigma = dt / (dx ** 2)
# предварительно составляем матрицу из нулей
nullMatrix = np.zeros((Nt, Nx))


print('Количество временных слоёв', Nt)
print('Шагов по пространству на одном временном слое', Nx)
print('Число Куранта', sigma)


# Создадим функцию, которая заполняет первую строку матрицы
def firstString():
    for j in range(0, Nx - 1):
        nullMatrix[0, j] = B1
    nullMatrix[0, Nx - 1] = B2


def Borders():
    for j in range(0, Nt):
        nullMatrix[j, 0] = B1
        nullMatrix[j, Nx - 1] = B2


# После заполнения первой строки, необходимо рассчитать прогоночные коэффициенты

def nextAlpha(preAlpha):
    result = Decimal(sigma / (1 + (2 * sigma) - (sigma * preAlpha)))
    return result


def nextBetta(preBetta, alpha, w):
    a = 8
    result = Decimal((alpha / sigma) * (w + sigma * preBetta - dt * (w * (1 - w) * (a - w))))
    return result


def ourFunc(alpha, betta, w):
    result = Decimal(alpha * w + betta)
    return result

def Method():
    for i in range(1, Nt - 1):
        Alphas = np.zeros(Nx)
        Alphas[1] = 0
        for j in range(1, Nx - 1):
            Alphas[j+1] = nextAlpha(Alphas[j])
        Bettas = np.zeros(Nx)
        Bettas[1] = B1
        for j in range(1, Nx - 1):
            Bettas[j + 1] = nextBetta(Bettas[j], Alphas[j + 1], nullMatrix[i - 1, j])
        for j in range(Nx - 2, 0, -1):
            nullMatrix[i + 1, j] = ourFunc(Alphas[j + 1], Bettas[j + 1], nullMatrix[i + 1, j + 1])
            
firstString()
Borders()

Method()
endTime = time.time() - origTime
print('Общее время расчёта =', endTime, 'секунд')
J = []
for j in range(0, len(nullMatrix[0])):
    J.append(j)
plt.ylabel('Особей, доля')
plt.xlabel('Ареал, Км')
plt.plot(J, nullMatrix[Nt - 1])
plt.show()
