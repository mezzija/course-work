import matplotlib.pyplot as plt
import matplotlib.animation as ani

t0=10
c1=2
c2=3
d=2
L = 5
c0 = L/2*d+(c1+c2)/2
k1=2
k2=3

def heat_capacity(t):
    if t < t0-d:
       return c1
    elif t > t0+d:
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
solve = []
x_max = 1
_time = 1

y0 = []
for i in range(0, x_nodes):
    y0.append(0)

solve.append(y0)

for i in range()