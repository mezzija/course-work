import matplotlib.pyplot as plt
import matplotlib.animation as ani


t0=10
c1=2
c2=20
d=0.01
L = 100
c0 = L/(2*d)+(c1+c2)/2
k1=1
k2=5

def heat_capacity(t):
    if t <= t0-d:
       return c1
    elif t >= t0+d:
        return c2
    elif t > t0-d and t < t0:
        return (t0 * c1 - (t0 - d) * c0 - (c1 - c0) * t) / (t0 - (t0 - d))
    elif t < t0+d and t>t0:
        return ((t0 + d) * c0 - t0 * c2 - (c0-c2)*t) / (t0 + d - t0)
    elif t == t0:
        return c0

def heat_conductivity(t):
    if t <= t0-d:
        return k1
    elif t > t0-d and t < t0+d:
        return ((t0 + d)*k1 - (t0-d)*k2-(k1-k2)*t)/(2*d)
    elif t >= t0+d:
        return k2

x_nodes = 100
t_nodes = 100
y_nodes = 100

y_max = 1
x_max = 1
_time = 1


left_x = 50
right_x = 50

left_y = 50
right_y = 50

dt = _time/t_nodes
dx = x_max/x_nodes
dy = y_max/y_nodes


solve = [0 for i in range(t_nodes)]

for i in range(t_nodes):
    solve[i] = [0 for j in range(x_nodes)]
    for j in range(x_nodes):
        solve[i][j] = [0 for j in range(y_nodes)]


# начальные условия
for i in range(x_nodes):
    for j in range(y_nodes):
        solve[0][i][j] = 0


for i in range(1, t_nodes):
    for j in range(y_nodes):
        alpha = [0 for k in range(x_nodes)]
        betta = [0 for k in range(x_nodes)]
        alpha[1] = 0
        betta[1] = left_x
        for k in range(2, x_nodes):
            y1 = (solve[i - 1][j][k-1]+solve[i-1][j][k-2])/2
            y2 = (solve[i - 1][j][k-1]+solve[i-1][j][k])/2
            a = -(dt*heat_conductivity(y2)/(heat_capacity(solve[i-1][j][k-1])*dx*dx))
            b = -(dt * heat_conductivity(y1) / (heat_capacity(solve[i - 1][j][k-1]) * dx * dx))
            c = (1 + (dt * heat_conductivity(y2) / (heat_capacity(solve[i - 1][j][k-1]) * dx * dx)) + (dt * heat_conductivity(y1) / (heat_capacity(solve[i - 1][j][k-1]) * dx * dx)))
            alpha[k] = -(b / (a * alpha[k - 1] + c))
            betta[k] = (solve[i - 1][j][k-1] - a * betta[k - 1]) / (a * alpha[k - 1] + c)

        solve[i][j][x_nodes-1] = right_x
        for k in range(x_nodes - 2, -1, -1):
            solve[i][j][k] = alpha[k + 1] * solve[i][j][k+1] + betta[k + 1]
    for j in range(x_nodes):
        alpha = [0 for k in range(y_nodes)]
        betta = [0 for k in range(y_nodes)]
        alpha[1] = 0
        betta[1] = left_y
        for k in range(2, y_nodes):
            y1 = (solve[i][k - 1][j] + solve[i][k - 2][j]) / 2
            y2 = (solve[i][k - 1][j] + solve[i][k][j]) / 2
            a = -(dt * heat_conductivity(y2) / (heat_capacity(solve[i][k - 1][j]) * dx * dx))
            b = -(dt * heat_conductivity(y1) / (heat_capacity(solve[i][k - 1][j]) * dx * dx))
            c = (1 + (dt * heat_conductivity(y2) / (heat_capacity(solve[i][k - 1][j]) * dx * dx)) + (dt * heat_conductivity(y1) / (heat_capacity(solve[i][k-1][j]) * dx * dx)))
            alpha[k] = -(b / (a * alpha[k - 1] + c))
            betta[k] = (solve[i][k-1][j] - a * betta[k - 1]) / (a * alpha[k - 1] + c)

        solve[i][y_nodes - 1][j] = right_y
        for k in range(y_nodes - 2, -1, -1):
            solve[i][k][j] = alpha[k + 1] * solve[i][k + 1][j] + betta[k + 1]




plt.imshow(solve[5])
plt.colorbar()
plt.show()
















