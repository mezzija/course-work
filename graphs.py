import matplotlib.pyplot as plt
import numpy as np
import time
import math

# количество узлов по времени
t_nodes = 100
# количество пространственных узлов
x_nodes = 20
y_nodes = 20
z_nodes = 20

with open('C:\\Users\mezzi\Desktop\курсач 2.0\course-work\curs.txt', 'r+') as f:
    a = np.loadtxt(f)
    solve = a.reshape((t_nodes, x_nodes, y_nodes, z_nodes))

with open('C:\\Users\mezzi\Desktop\курсач 2.0\course-work\cca.txt', 'r+') as f:
    array = np.loadtxt(f)
    exact_solution = array.reshape((t_nodes, x_nodes, y_nodes, z_nodes))

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
plt.show()
