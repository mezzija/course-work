import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np
import time
import math

# температура плавления
t0 = 0
# температура в центре
U0 = 1

k1 = 0.5
k2 = 0.75
c1 = 2
c2 = 1.25
# энтальпия
L = 1
# дельта
d = 0.15

c0 = L / (2 * d) + (c1 + c2) / 2


def heat_capacity(t):
    if t <= t0 - d:
        return c1
    elif t >= t0 + d:
        return c2
    elif t > t0 - d and t < t0:
        return (t0 * c1 - (t0 - d) * c0 - (c1 - c0) * t) / (t0 - (t0 - d))
    elif t < t0 + d and t > t0:
        return ((t0 + d) * c0 - t0 * c2 - (c0 - c2) * t) / (t0 + d - t0)
    elif t == t0:
        return c0


def heat_conductivity(t):
    if t <= t0 - d:
        return k1
    elif t > t0 - d and t < t0 + d:
        return ((t0 + d) * k1 - (t0 - d) * k2 - (k1 - k2) * t) / (2 * d)
    elif t >= t0 + d:
        return k2


# количество узлов по времени
t_nodes = 10
# количество пространственных узлов
x_nodes = 10
y_nodes = 10
z_nodes = 10
# максимальная длинна длинна по координатам
y_max = 4
x_max = 4
z_max = 4
# максимальное время
t_max = 64
# шаги
dt = t_max / t_nodes
dx = x_max / x_nodes
dy = y_max / y_nodes
dz = z_max / z_nodes
# средние значения
x0 = x_max / 2
y0 = y_max / 2
z0 = z_max / 2

# коэфициент при кси
a = 0.2


# радиус жидкой фазы
def ksi(t):
    return a * math.sqrt(t_max - t)


# расчет радиуса
def r2(x, y, z):
    return (x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2


# коэфциенты для точного решения
A1 = (k2 / k1) * ((U0 - t0) / (a ** 2)) + (a * L) / (4 * k1)
A2 = (U0 - t0) / a ** 2
B1 = t0 + a ** 2 * A1
B2 = U0


# рассчет температуры
def U(t, x, y, z):
    if math.sqrt(r2(x, y, z)) >= ksi(t):
        return B1 - A1 * (r2(x, y, z) / (t0 - t))
    elif math.sqrt(r2(x, y, z)) < ksi(t):
        return B2 - A2 * (r2(x, y, z) / (t0 - t))
