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
x_nodes = 100
y_nodes = 100
z_nodes = 100
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


# расчет температуры
def U(t, x, y, z):
    if math.sqrt(r2(x, y, z)) >= ksi(t):
        return B1 - A1 * (r2(x, y, z) / (t_max - t))
    elif math.sqrt(r2(x, y, z)) < ksi(t):
        return B2 - A2 * (r2(x, y, z) / (t_max - t))


# расчет функции
def f(t, x, y, z):
    if math.sqrt(r2(x, y, z)) >= ksi(t):
        return A1 * ((6 * k1) / (t_max - t) - (c1 * r2(x, y, z)) / (t_max - t) ** 2)
    elif math.sqrt(r2(x, y, z)) < ksi(t):
        return A2 * ((6 * k2) / (t_max - t) - (c2 * r2(x, y, z)) / (t_max - t) ** 2)


# список решения
solve = []
# начальные условия
initial_conditions = np.zeros(shape=(x_nodes, y_nodes, z_nodes), dtype=float)
for x in range(x_nodes):
    for y in range(y_nodes):
        for z in range(z_nodes):
            initial_conditions[x, y, z] = U(0, x * dx, y * dy, z * dz)

solve.append(initial_conditions)


start_time = time.time()
for t in range(1, t_nodes):
    interval_matrix_X = np.zeros(shape=(z_nodes, y_nodes, x_nodes), dtype=float)
    interval_matrix_Y = np.zeros(shape=(z_nodes, y_nodes, x_nodes), dtype=float)
    layer_matrix = np.zeros(shape=(z_nodes, y_nodes, x_nodes), dtype=float)
    for y in range(y_nodes):
        for z in range(z_nodes):
            alpha = [0] * x_nodes
            betta = [0] * x_nodes
            alpha[1] = 0
            betta[1] = U(t, 0, y * dy, z * dz)
            for x in range(2, x_nodes):
                y1 = (solve[t - 1][x - 1][y][z] + solve[t - 1][x - 2][y][z]) / 2
                y2 = (solve[t - 1][x - 1][y][z] + solve[t - 1][x][y][z]) / 2
                a = -(dt * heat_conductivity(y2) / (heat_capacity(solve[t - 1][x - 1][y][z]) * dx * dx))
                b = -(dt * heat_conductivity(y1) / (heat_capacity(solve[t - 1][x - 1][y][z]) * dx * dx))
                c = (1 + (dt * heat_conductivity(y2) / (heat_capacity(solve[t - 1][x - 1][y][z]) * dx * dx)) + (
                        dt * heat_conductivity(y1) / (heat_capacity(solve[t - 1][x - 1][y][z]) * dx * dx)))
                alpha[x] = -(b / (a * alpha[x - 1] + c))
                betta[x] = (solve[t - 1][x - 1][y][z] - a * betta[x - 1]) / (a * alpha[x - 1] + c)
            interval_matrix_X[x_nodes - 1][y][z] = U(t, x_max, y * dy, z * dz)
            for x in range(x_nodes - 2, -1, -1):
                interval_matrix_X[x][y][z] = alpha[x + 1] * interval_matrix_X[x + 1][y][z] + betta[x + 1]
    for x in range(x_nodes):
        for z in range(z_nodes):
            alpha = [0 for a in range(y_nodes)]
            betta = [0 for b in range(y_nodes)]
            alpha[1] = 0
            betta[1] = U(t, x * dx, 0, z * dz)
            for y in range(2, y_nodes):
                y1 = (interval_matrix_X[x][y - 1][z] + interval_matrix_X[x][y - 2][z]) / 2
                y2 = (interval_matrix_X[x][y - 1][z] + interval_matrix_X[x][y][z]) / 2
                a = -(dt * heat_conductivity(y2) / (heat_capacity(interval_matrix_X[x][y - 1][z]) * dx * dx))
                b = -(dt * heat_conductivity(y1) / (heat_capacity(interval_matrix_X[x][y - 1][z]) * dx * dx))
                c = (1 + (dt * heat_conductivity(y2) / (heat_capacity(interval_matrix_X[x][y - 1][z]) * dx * dx)) + (
                        dt * heat_conductivity(y1) / (heat_capacity(interval_matrix_X[x][y - 1][z]) * dx * dx)))
                alpha[y] = -(b / (a * alpha[y - 1] + c))
                betta[y] = (interval_matrix_X[x][y - 1][z] - a * betta[y - 1]) / (a * alpha[y - 1] + c)
            interval_matrix_Y[x][y_nodes - 1][z] = U(t, x * dx, y_max, z * dz)
            for y in range(y_nodes - 2, -1, -1):
                interval_matrix_Y[x][y][z] = alpha[y + 1] * interval_matrix_Y[x][y + 1][z] + betta[y + 1]
    for x in range(x_nodes):
        for y in range(y_nodes):
            alpha = [0 for a in range(z_nodes)]
            betta = [0 for b in range(z_nodes)]
            alpha[1] = 0
            betta[1] = U(t, x * dx, y * dy, 0)
            for z in range(2, z_nodes):
                y1 = (interval_matrix_Y[x][y][z - 1] + interval_matrix_Y[x][y][z - 2]) / 2
                y2 = (interval_matrix_Y[x][y][z - 1] + interval_matrix_Y[x][y][z]) / 2
                a = -(dt * heat_conductivity(y2) / (heat_capacity(interval_matrix_Y[x][y][z - 1]) * dx * dx))
                b = -(dt * heat_conductivity(y1) / (heat_capacity(interval_matrix_Y[x][y][z - 1]) * dx * dx))
                c = (1 + (dt * heat_conductivity(y2) / (heat_capacity(interval_matrix_Y[x][y][z - 1]) * dx * dx)) + (
                        dt * heat_conductivity(y1) / (heat_capacity(interval_matrix_Y[x][y][z - 1]) * dx * dx)))
                alpha[z] = -(b / (a * alpha[z - 1] + c))
                betta[z] = (interval_matrix_Y[x][y][z - 1] - a * betta[z - 1]) / (a * alpha[z - 1] + c)
            layer_matrix[x][y][z_nodes - 1] = U(t, x * dx, y * dy, z_max)
            for z in range(z_nodes - 2, -1, -1):
                layer_matrix[x][y][z] = alpha[z + 1] * interval_matrix_Y[x][y][z + 1] + betta[z + 1]
    solve.append(layer_matrix)


plt.imshow(solve[1][5])
plt.colorbar()
print("--- %s seconds ---" % (time.time() - start_time))
plt.show()

