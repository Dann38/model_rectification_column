import numpy as np
from scipy.integrate import quad
from sol import sol
import matplotlib.pyplot as plt

s = [0, 2]
t = [0, 4]
c = [1, 3]
m = 29
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


s = np.linspace(0, 2, 29)
x = s**2+np.cos(t[1])
y = np.exp(s)*(np.sin(t[1])+2*np.cos(t[1]))

fig, axs = plt.subplots(nrows= 2 , ncols= 1 )
axs[0].set_title("Численное решение (узлов по t: 232, узлов по s: 29)")
axs[0].plot(s, x)
axs[0].plot(rez["s"], rez["x(t1)"], "o")
axs[0].grid()
axs[0].legend(["аналитическое решение", "численное решение"])
axs[0].set_ylabel("$x(s, t_1)$")

axs[1].plot(s, y)
axs[1].plot(rez["s"], rez["y(t1)"], "o")
axs[1].grid()
axs[1].legend(["аналитическое решение", "численное решение"])
axs[1].set_ylabel("$y(s, t_1)$")
axs[1].set_xlabel("$s$")
plt.show()
print(f"max x:{np.max(abs(rez['s']**2+np.cos(t[1])-rez['x(t1)']))}")
print(f"max y:{np.max(abs(np.exp(rez['s'])*(np.sin(t[1])+2*np.cos(t[1]))-rez['y(t1)']))}")
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
