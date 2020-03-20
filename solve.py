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
solve = [0 for i in range(t_nodes)]
x_max = 1
_time = 1
left = 0
right = 20
dt = _time/t_nodes
dx = x_max/x_nodes

for i in range(t_nodes):
    solve[i] = [0 for j in range(t_nodes)]

for i in range(x_nodes):
    solve[0][i] = 0

solve[0][0] = 0
solve[0][x_nodes-1] = 20


for i in range(1, t_nodes-1):
    alpha = [0 for k in range(100)]
    betta = [0 for k in range(100)]
    alpha[1] = 0
    betta[1] = left
    for j in range(2, x_nodes):
        y1 = (solve[i-1][j-1]+solve[i-1][j-2])/2
        y2 = (solve[i-1][j-1]+solve[i-1][j])/2
        a = -(dt*heat_conductivity(y2)/(heat_capacity(solve[i-1][j-1])*dx*dx))
        b = -(dt*heat_conductivity(y1)/(heat_capacity(solve[i-1][j-1])*dx*dx))
        c = (1+(dt*heat_conductivity(y2)/(heat_capacity(solve[i-1][j-1])*dx*dx))+(dt*heat_conductivity(y1)/(heat_capacity(solve[i-1][j-1])*dx*dx)))
        alpha[j] = -(b/(a*alpha[j-1]+c))
        betta[j] = (solve[i-1][j]-a*betta[j-1])/(a*alpha[j-1]+c)

    solve[i][x_nodes-1] = right
    for j in range(x_nodes-2, 0, -1):
        solve[i][j] = alpha[j+1]*solve[i][j+1]+betta[j+1]

fig = plt.figure()
def animate(n):
    line = plt.plot(solve[n], color='g')
    return line

a = ani.FuncAnimation(fig, animate, frames=len(solve), interval=100, blit=True, repeat=True)
plt.show()
