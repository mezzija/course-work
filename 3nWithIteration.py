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
t_nodes = 100
# количество пространственных узлов
x_nodes = 20
y_nodes = 20
z_nodes = 20
# максимальная длинна длинна по координатам
y_max = 4
x_max = 4
z_max = 4
# максимальное время
t_max = 64
# шаги
dt = 0.5
dx = x_max / x_nodes
dy = y_max / y_nodes
dz = z_max / z_nodes
# средние значения
x0 = x_max / 2
y0 = y_max / 2
z0 = z_max / 2

# коэфициент при кси
aa = 0.2


# радиус жидкой фазы
def ksi(t):
    return aa * math.sqrt(t_max - t)


# расчет радиуса
def r2(x, y, z):
    return ((x - x0) ** 2) + ((y - y0) ** 2) + ((z - z0) ** 2)


# коэфциенты для точного решения
A1 = (k2 / k1) * ((U0 - t0) / (aa ** 2)) + (aa * L) / (4 * k1)
A2 = (U0 - t0) / (aa ** 2)
B1 = t0 + (aa ** 2) * A1
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
        return A1 * ((6 * k1) / (t_max - t) - (c1 * r2(x, y, z)) / ((t_max - t) ** 2))
    elif math.sqrt(r2(x, y, z)) < ksi(t):
        return A2 * ((6 * k2) / (t_max - t) - (c2 * r2(x, y, z)) / ((t_max - t) ** 2))


# список решения
solve = np.zeros(shape=(t_nodes, x_nodes, y_nodes, z_nodes), dtype=float)
# начальные условия
initial_conditions = np.zeros(shape=(x_nodes, y_nodes, z_nodes), dtype=float)
for x in range(x_nodes):
    for y in range(y_nodes):
        for z in range(z_nodes):
            initial_conditions[x, y, z] = U(0, x * dx, y * dy, z * dz)

# solve.append(initial_conditions)
solve[0] = initial_conditions

exact_solution = np.zeros(shape=(t_nodes, x_nodes, y_nodes, z_nodes), dtype=float)
for t in range(t_nodes):
    for x in range(x_nodes):
        for y in range(y_nodes):
            for z in range(z_nodes):
                exact_solution[t, x, y, z] = U(t * dt, x * dx, y * dy, z * dz)

start_time = time.time()
for t in range(1, t_nodes):
    interval_matrix_X = np.zeros(shape=(x_nodes, y_nodes, z_nodes), dtype=float)
    interval_matrix_Y = np.zeros(shape=(x_nodes, y_nodes, z_nodes), dtype=float)
    layer_matrix = np.zeros(shape=(x_nodes, y_nodes, z_nodes), dtype=float)
    iteration_matrix = np.zeros(shape=(x_nodes, y_nodes, z_nodes), dtype=float)
    for x in range(x_nodes):
        for y in range(y_nodes):
            for z in range(z_nodes):
                iteration_matrix[x][y][z] = solve[t - 1][x][y][z]
    n = 0
    while True:
        for y in range(y_nodes):
            for z in range(z_nodes):
                alpha = [0 for a in range(x_nodes)]
                betta = [0 for a in range(x_nodes)]
                alpha[1] = 0
                betta[1] = U(t * dt, 0, y * dy, z * dz)
                #interval_matrix_X[0][y][z] = U(t * dt, 0, y * dy, z * dz)
                for x in range(2, x_nodes):
                    y1 = (iteration_matrix[x - 1][y][z] + iteration_matrix[x - 2][y][z]) / 2
                    y2 = (iteration_matrix[x - 1][y][z] + iteration_matrix[x][y][z]) / 2
                    A = -(dt * heat_conductivity(y2) / (heat_capacity(iteration_matrix[x - 1][y][z]) * dx * dx))
                    B = -(dt * heat_conductivity(y1) / (heat_capacity(iteration_matrix[x - 1][y][z]) * dx * dx))
                    C = (1 + (dt * heat_conductivity(y2) / (heat_capacity(iteration_matrix[x - 1][y][z]) * dx * dx)) + (
                            dt * heat_conductivity(y1) / (heat_capacity(iteration_matrix[x - 1][y][z]) * dx * dx)))
                    F = iteration_matrix[x - 1][y][z] + ((f(t * dt, (x - 1) * dx, y * dy, z * dz) / 3)*dt)/heat_capacity(iteration_matrix[x - 1][y][z])
                    alpha[x] = -B / (C + A * alpha[x - 1])
                    betta[x] = (F - A * betta[x - 1]) / (C + A * alpha[x - 1])
                interval_matrix_X[x_nodes - 1][y][z] = U(t * dt, x_max, y * dy, z * dz)
                for x in range(x_nodes - 2, -1, -1):
                    interval_matrix_X[x][y][z] = alpha[x + 1] * interval_matrix_X[x + 1][y][z] + betta[x + 1]
        for x in range(x_nodes):
            for z in range(z_nodes):
                alpha = [0 for a in range(y_nodes)]
                betta = [0 for a in range(y_nodes)]
                alpha[1] = 0
                betta[1] = U(t * dt, x * dx, 0, z * dz)
                #interval_matrix_Y[x][0][z] = U(t * dt, x * dx, 0, z * dz)
                for y in range(2, y_nodes):
                    y1 = (interval_matrix_X[x][y - 1][z] + interval_matrix_X[x][y - 2][z]) / 2
                    y2 = (interval_matrix_X[x][y - 1][z] + interval_matrix_X[x][y][z]) / 2
                    A = -(dt * heat_conductivity(y2) / (heat_capacity(interval_matrix_X[x][y - 1][z]) * dx * dx))
                    B = -(dt * heat_conductivity(y1) / (heat_capacity(interval_matrix_X[x][y - 1][z]) * dx * dx))
                    C = (1 + (dt * heat_conductivity(y2) / (
                            heat_capacity(interval_matrix_X[x][y - 1][z]) * dx * dx)) + (
                                 dt * heat_conductivity(y1) / (
                                 heat_capacity(interval_matrix_X[x][y - 1][z]) * dx * dx)))
                    F = interval_matrix_X[x][y - 1][z] + ((f(t * dt, x * dx, (y - 1) * dy, z * dz) / 3)*dt)/heat_capacity(interval_matrix_X[x][y - 1][z])
                    alpha[y] = -B / (C + A * alpha[y - 1])
                    betta[y] = (F - A * betta[y - 1]) / (C + A * alpha[y - 1])
                interval_matrix_Y[x][y_nodes - 1][z] = U(t * dt, x * dx, y_max, z * dz)
                for y in range(y_nodes - 2, -1, -1):
                    interval_matrix_Y[x][y][z] = alpha[y + 1] * interval_matrix_Y[x][y + 1][z] + betta[y + 1]
        for x in range(x_nodes):
            for y in range(y_nodes):
                alpha = [0 for a in range(z_nodes)]
                betta = [0 for b in range(z_nodes)]
                alpha[1] = 0
                betta[1] = U(t * dt, x * dx, y * dy, 0)
                #layer_matrix[x][y][0] = U(t * dt, x * dx, y * dy, 0)
                for z in range(2, z_nodes):
                    y1 = (interval_matrix_Y[x][y][z - 1] + interval_matrix_Y[x][y][z - 2]) / 2
                    y2 = (interval_matrix_Y[x][y][z - 1] + interval_matrix_Y[x][y][z]) / 2
                    A = -(dt * heat_conductivity(y2) / (heat_capacity(interval_matrix_Y[x][y][z - 1]) * dx * dx))
                    B = -(dt * heat_conductivity(y1) / (heat_capacity(interval_matrix_Y[x][y][z - 1]) * dx * dx))
                    C = (1 + (dt * heat_conductivity(y2) / (
                            heat_capacity(interval_matrix_Y[x][y][z - 1]) * dx * dx)) + (
                                 dt * heat_conductivity(y1) / (
                                 heat_capacity(interval_matrix_Y[x][y][z - 1]) * dz ** 2)))
                    F = interval_matrix_Y[x][y][z - 1] + ((f(t * dt, x * dx, y * dy, (z - 1) * dz) / 3)*dt)/heat_capacity(interval_matrix_Y[x][y][z - 1])
                    alpha[z] = -B / (C + A * alpha[z - 1])
                    betta[z] = (F - A * betta[z - 1]) / (C + A * alpha[z - 1])
                layer_matrix[x][y][z_nodes - 1] = U(t * dt, x * dx, y * dy, z_max)
                for z in range(z_nodes - 2, -1, -1):
                    layer_matrix[x][y][z] = alpha[z + 1] * interval_matrix_Y[x][y][z + 1] + betta[z + 1]
        mistake = np.zeros(shape=(x_nodes, y_nodes, z_nodes), dtype=float)

        for x in range(x_nodes):
            for y in range(y_nodes):
                for z in range(z_nodes):
                    a = math.fabs(iteration_matrix[x][y][z])
                    b = math.fabs(layer_matrix[x, y, z])
                    mistake[x, y, z] = float(math.fabs(a - b))

        err = mistake.max()
        print(err)
        if n > 15:
            break
        elif err > 0.001:
            for x in range(x_nodes):
                for y in range(y_nodes):
                    for z in range(z_nodes):
                        iteration_matrix[x, y, z] = layer_matrix[x, y, z]
            n = n + 1
        elif err < 0.001:
            print(n)
            break
    solve[t] = layer_matrix

'''
array = solve.reshape((t_nodes, x_nodes, y_nodes, z_nodes))
with open('C:\\Users\mezzi\Desktop\курсач 2.0\course-work\curs.txt', 'w+') as f:
    for each_3d in array:
        for each_2d in each_3d:
            np.savetxt(f, each_2d)

array = exact_solution.reshape((t_nodes, x_nodes, y_nodes, z_nodes))
with open('C:\\Users\mezzi\Desktop\курсач 2.0\course-work\cca.txt', 'w+') as f:
    for each_3d in array:
        for each_2d in each_3d:
            np.savetxt(f, each_2d)

print("--- %s seconds ---" % (time.time() - start_time))

'''
inconsistency = np.zeros(shape=(t_nodes, x_nodes, y_nodes, z_nodes), dtype=float)
for t in range(t_nodes):
    for x in range(x_nodes):
        for y in range(y_nodes):
            for z in range(z_nodes):
                inconsistency[t, x, y, z] = solve[t][x][y][z] - exact_solution[t][x][y][z]

mistake = []
for t in range(t_nodes):
    interval = np.zeros(shape=(x_nodes, y_nodes, z_nodes))
    for x in range(x_nodes):
        for y in range(y_nodes):
            for z in range(z_nodes):
                interval[x, y, z] = math.fabs(inconsistency[t, x, y, z])
    mistake.append(interval.max())
print(mistake)
plt.plot(mistake)
plt.show()

plt.imshow(inconsistency[5][5])
plt.colorbar()
plt.show()

plt.imshow(exact_solution[5][5])
plt.colorbar()
plt.show()

plt.imshow(solve[5][5])
plt.colorbar()
print("--- %s seconds ---" % (time.time() - start_time))
plt.show()
