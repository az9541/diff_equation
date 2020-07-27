from math import floor
import numpy as np
import matplotlib.pyplot as plt
import time
from decimal import Decimal

# изначально flops = 0
flops = 0
# начинаем записывать время
origTime = time.time()
# задаём граничные условия
B1 = np.double(10)  # B(Borders)
B2 = np.double(0)
# шаги по пространству и по времени
dx = np.double(0.1)
dt = np.double(0.000005)
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

print('Количество временных слоёв:', Nt)
print('Шагов по пространству на одном временном слое:', Nx)
print('Число Куранта:', sigma)


# заполняем первую строчку матрицы левыми граничными условиями
def firstString():
    for j in range(0, Nx - 1):
        nullMatrix[0, j] = B1


# задаём правые и левые граничные условия на всём промежутке времени
def Borders():
    for j in range(0, Nt):
        nullMatrix[j, 0] = B1
        nullMatrix[j, Nx - 1] = B2


# записываем разностную схему уравнения в явном виде. a - какая-то константа, FLOPS = 22

def Wi(i, j, w):
    a = 4
    solution = sigma * ((w[i, j + 1]) - 2 * w[i, j] + w[i, j - 1]) + w[i, j] - dt * w[i, j] * (1 - w[i, j]) * (a - w[i, j])
    return Decimal(solution)


# заполняем первую строку и применяем ГУ
firstString()
Borders()
# для каждого временного для каждого элемента строки применяем функцию Wi
for i in range(1, Nt):
    for j in range(1, Nx - 1):
        L = Wi(i - 1, j, nullMatrix)
        nullMatrix[i, j] = L
        flops += 22
# создаём вспомогательный массив для построения графика
J = []
numberOfDotsX = len(nullMatrix[Nt - 1])
for j in range(0, numberOfDotsX):
    J.append(j)
# заканчиваем записывать время расчёта
endTime = time.time() - origTime
print('Общее время расчёта =', endTime, 'секунд')
print('Количество FLOPS =', flops)
# строим график по последней строке
plt.plot(J, nullMatrix[Nt - 1], color = 'red')
plt.show()