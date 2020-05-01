import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np
import time

t0 = 10
c1 = 2
c2 = 20
d = 0.01
L = 100
c0 = L / (2 * d) + (c1 + c2) / 2
k1 = 1
k2 = 5


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


t_nodes = 10

x_nodes = 10
y_nodes = 10
z_nodes = 10

y_max = 1
x_max = 1
z_max = 1

_time = 1

left_x = 50
right_x = 50

left_y = 50
right_y = 50

left_z = 50
right_z = 50

dt = _time / t_nodes
dx = x_max / x_nodes
dy = y_max / y_nodes
dz = z_max / z_nodes

solve = []
initial_conditions = np.zeros(shape=(z_nodes, y_nodes, x_nodes), dtype=float)
initial_conditions.fill(1)
solve.append(initial_conditions)

start_time = time.time()
for t in range(1, t_nodes):
    interval_matrix_X = np.zeros(shape=(z_nodes, y_nodes, x_nodes), dtype=float)
    interval_matrix_X.fill(0)
    interval_matrix_Y = np.zeros(shape=(z_nodes, y_nodes, x_nodes), dtype=float)
    interval_matrix_Y.fill(0)
    layer_matrix = np.zeros(shape=(z_nodes, y_nodes, x_nodes), dtype=float)
    layer_matrix.fill(0)

    for y in range(y_nodes):
        for z in range(z_nodes):
            alpha = [0 for a in range(x_nodes)]
            betta = [0 for b in range(x_nodes)]
            alpha[1] = 0
            betta[1] = left_x
            for x in range(2, x_nodes):
                y1 = (solve[t - 1][x - 1][y][z] + solve[t - 1][x - 2][y][z]) / 2
                y2 = (solve[t - 1][x - 1][y][z] + solve[t - 1][x][y][z]) / 2
                a = -(dt * heat_conductivity(y2) / (heat_capacity(solve[t - 1][x - 1][y][z]) * dx * dx))
                b = -(dt * heat_conductivity(y1) / (heat_capacity(solve[t - 1][x - 1][y][z]) * dx * dx))
                c = (1 + (dt * heat_conductivity(y2) / (heat_capacity(solve[t - 1][x - 1][y][z]) * dx * dx)) + (
                        dt * heat_conductivity(y1) / (heat_capacity(solve[t - 1][x - 1][y][z]) * dx * dx)))
                alpha[x] = -(b / (a * alpha[x - 1] + c))
                betta[x] = (solve[t - 1][x - 1][y][z] - a * betta[x - 1]) / (a * alpha[x - 1] + c)
            interval_matrix_X[x_nodes - 1][y][z] = right_x

            for x in range(x_nodes - 2, -1, -1):
                interval_matrix_X[x][y][z] = alpha[x + 1] * interval_matrix_X[x + 1][y][z] + betta[x + 1]
    for x in range(x_nodes):
        for z in range(z_nodes):
            alpha = [0 for a in range(y_nodes)]
            betta = [0 for b in range(y_nodes)]
            alpha[1] = 0
            betta[1] = left_y
            for y in range(2, y_nodes):
                y1 = (interval_matrix_X[x][y - 1][z] + interval_matrix_X[x][y - 2][z]) / 2
                y2 = (interval_matrix_X[x][y - 1][z] + interval_matrix_X[x][y][z]) / 2
                a = -(dt * heat_conductivity(y2) / (heat_capacity(interval_matrix_X[x][y - 1][z]) * dx * dx))
                b = -(dt * heat_conductivity(y1) / (heat_capacity(interval_matrix_X[x][y - 1][z]) * dx * dx))
                c = (1 + (dt * heat_conductivity(y2) / (heat_capacity(interval_matrix_X[x][y - 1][z]) * dx * dx)) + (
                        dt * heat_conductivity(y1) / (heat_capacity(interval_matrix_X[x][y - 1][z]) * dx * dx)))
                alpha[y] = -(b / (a * alpha[y - 1] + c))
                betta[y] = (interval_matrix_X[x][y - 1][z] - a * betta[y - 1]) / (a * alpha[y - 1] + c)
            interval_matrix_Y[x][y_nodes - 1][z] = right_y
            for y in range(y_nodes - 2, -1, -1):
                interval_matrix_Y[x][y][z] = alpha[y + 1] * interval_matrix_Y[x][y + 1][z] + betta[y + 1]

    for x in range(x_nodes):
        for y in range(y_nodes):
            alpha = [0 for a in range(z_nodes)]
            betta = [0 for b in range(z_nodes)]
            alpha[1] = 0
            betta[1] = left_z
            for z in range(2, z_nodes):
                y1 = (interval_matrix_Y[x][y][z - 1] + interval_matrix_Y[x][y][z - 2]) / 2
                y2 = (interval_matrix_Y[x][y][z - 1] + interval_matrix_Y[x][y][z]) / 2
                a = -(dt * heat_conductivity(y2) / (heat_capacity(interval_matrix_Y[x][y][z - 1]) * dx * dx))
                b = -(dt * heat_conductivity(y1) / (heat_capacity(interval_matrix_Y[x][y][z - 1]) * dx * dx))
                c = (1 + (dt * heat_conductivity(y2) / (heat_capacity(interval_matrix_Y[x][y][z - 1]) * dx * dx)) + (
                        dt * heat_conductivity(y1) / (heat_capacity(interval_matrix_Y[x][y][z - 1]) * dx * dx)))
                alpha[z] = -(b / (a * alpha[z - 1] + c))
                betta[z] = (interval_matrix_Y[x][y][z - 1] - a * betta[z - 1]) / (a * alpha[z - 1] + c)
            layer_matrix[x][y][z_nodes - 1] = right_z
            for z in range(z_nodes - 2, -1, -1):
                layer_matrix[x][y][z] = alpha[z + 1] * interval_matrix_Y[x][y][z + 1] + betta[z + 1]

    solve.append(layer_matrix)

plt.imshow(solve[1][5])
plt.colorbar()
print("--- %s seconds ---" % (time.time() - start_time))
plt.show()
