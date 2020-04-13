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

t_nodes = 100

x_nodes = 100
y_nodes = 100
z_nodes = 100

y_max = 1
x_max = 1
z_max = 1

_time = 1


left_x = 50
right_x = 50

left_y = 50
right_y = 50

left_z = 50
right_z =50

dt = _time/t_nodes
dx = x_max/x_nodes
dy = y_max/y_nodes
dz = z_max/z_nodes

solve = []


