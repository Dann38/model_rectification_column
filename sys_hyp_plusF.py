import numpy as np
from scipy.integrate import quad
from sol import sol
import matplotlib.pyplot as plt

s = [0, 2]
t = [0, 4]
c = [1, 3]
m = 30
B11 = lambda s, t: 1
B12 = lambda s, t: - 1/np.exp(s)
B21 = lambda s, t: np.exp(s)
B22 = lambda s, t: c[1]
B = [[B11, B12],
     [B21, B22]]

G11 = lambda t:  np.sin(t)/(np.exp(s[1])*(np.sin(t)+2*np.cos(t))-s[1]**2-np.cos(t))
G12 = lambda t: - G11(t)
G21 = lambda t: (2*np.sin(t) - np.cos(t))/(np.sin(t)+np.cos(t))
G22 = lambda t: - G21(t)

G = [[G11, G12],
     [G21, G22]]

F1 = lambda s, t: -2*c[0]*s-s**2+np.cos(t)
F2 = lambda s, t: -2*np.exp(s)*np.sin(t)-np.exp(s)*s**2
F = [F1, F2]

x0 = lambda s: s**2+1
y0 = lambda s: 2*np.exp(s)
fun0 = [x0, y0]

# s, c, fun0, B, F, G, m, Time_end
rez = sol(s, c, fun0, B, F, G, m, t[1])


s = np.linspace(0, 2, 30)
x = s**2+np.cos(t[1])
y = np.exp(s)*(np.sin(t[1])+2*np.cos(t[1]))
plt.subplot()
plt.plot(s, x)
plt.plot(rez["s"], rez["x(t1)"])
plt.ylim([-1, 4])
plt.subplots()
plt.plot(s, y)
plt.plot(rez["s"], rez["y(t1)"])
plt.ylim([-20, 0])
plt.show()

# t = np.linspace(t[0], t[-1], 60)
# xl = s[0]**2+np.cos(t)
# yl = np.exp(s[0])*(np.sin(t)+2*np.cos(t))
#
# plt.subplots()
# plt.plot(t, xl)
# plt.plot(rez["t(s0)"], rez["x(s0, t)"])
#
# plt.subplots()
# plt.plot(t, yl)
# plt.plot(rez["t(s0)"], rez["y(s0, t)"])
#
# plt.show()
