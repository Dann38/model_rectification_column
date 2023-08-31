import numpy as np
from scipy.integrate import quad
from sol import sol
import matplotlib.pyplot as plt

s = [0, 2]
t = [0, 4]
c = [1, 3]
m = 29
B11 = lambda s, t: -c[0]
B12 = lambda s, t: - np.exp(s)/(s+2)
B21 = lambda s, t: (s+2)/np.exp(s)
B22 = lambda s, t: c[1]/(s+2)

B = [[B11, B12],
     [B21, B22]]

G11 = lambda t:  np.exp(s[1])*np.sin(t)/((s[1]+2)*np.sin(t)-np.exp(s[1])*np.cos(t))
G12 = lambda t: - G11(t)
G21 = lambda t: 2*np.cos(t)/(np.cos(t)-2*np.sin(t))
G22 = lambda t: - G21(t)

G = [[G11, G12],
     [G21, G22]]

F1 = lambda s, t: 0
F2 = lambda s, t: 0
F = [F1, F2]

x0 = lambda s: np.exp(s)
y0 = lambda s: 0
fun0 = [x0, y0]



# s, c, fun0, B, F, G, m, Time_end
rez = sol(s, c, fun0, B, F, G, m, t[1])


s = np.linspace(0, 2, 30)
x = np.exp(s)*np.cos(t[1])
y = (s+2)*np.sin(t[1])
plt.subplot()
plt.plot(s, x)
plt.plot(rez["s"], rez["x(t1)"])

plt.subplots()
plt.plot(s, y)
plt.plot(rez["s"], rez["y(t1)"])

print(f"max x:{np.max(abs(np.exp(rez['s'])*np.cos(t[1])-rez['x(t1)']))}")
print(f"max y:{np.max(abs((rez['s']+2)*(np.sin(t[1]))-rez['y(t1)']))}")

plt.show()