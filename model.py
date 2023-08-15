import numpy as np
from scipy.integrate import quad
from sol import sol
import matplotlib.pyplot as plt

def related_task(s_, c, fun0, B, F, G, m, Time_end, phi):
    # Задача решается в обратном времени
    B11_inv = lambda s, t: -B[0][0](s, Time_end - t)
    B12_inv = lambda s, t: -B[1][0](s, Time_end - t)
    B21_inv = lambda s, t: -B[0][1](s, Time_end - t)
    B22_inv = lambda s, t: -B[1][1](s, Time_end - t)

    G11_inv = lambda t: 0
    G12_inv = lambda t: 0
    G21_inv = lambda t: 0
    G22_inv = lambda t: 0

    B_inv = [[B11_inv, B12_inv],
             [B21_inv, B22_inv]]
    G_inv = [[G11_inv, G12_inv],
             [G21_inv, G22_inv]]
    F_inv = [
        lambda s, t: F[0](s, Time_end - t),
        lambda s, t: F[1](s, Time_end - t)
    ]
    phi_0 = [
        lambda s: 0,
        lambda s: 0
    ]
    rez = sol(s_, c, phi_0, B_inv, F_inv, G_inv, m, Time_end)
    rez["x(t1)"] = list(reversed(rez["x(t1)"]))
    rez["y(t1)"] = list(reversed(rez["y(t1)"]))
    return rez

t1 = 4
t0 = 0
s0 = 0
s1 = 2
m = 200

c1 = 1
c2 = 3
A1 = lambda s, t: c1
A2 = lambda s, t: -np.exp(s)*(t1-t)/(s+2)
B1 = lambda s, t: 1/np.exp(s)
B2 = lambda s, t: -c2/(s+2)

Fx = lambda s, t: - np.exp(s)*np.cos(t1-t)
Fy = lambda s, t: (t1 - t)*np.cos(t1 - t) - (s+2)*np.cos(t1 - t)

G = None # TODO Подставить G
fun0 = None # TODO Подставить fun_0

s = [s0, s1]
c = [c1, c2]
B = [[A1, B1],
     [A2, B2]]
F = [Fx, Fy]

rez = related_task(s, c, fun0, B, F, G, m, t1, None)



s = np.linspace(0, 2, 30)
x = np.exp(s)*(t1-t0)*np.cos(t1-t0)
y = (s+2)*np.sin(t1-t0)
plt.subplot()
plt.plot(s, x)
plt.plot(rez["s"], rez["x(t1)"])
plt.subplots()
plt.plot(s, y)
plt.plot(rez["s"], rez["y(t1)"])
plt.show()

















# # Параметры =====================================
# c1 = 36  # м/ч L(s, t)/Hx(s, t)
# c2 = 520  # м/ч V(s, t)/Hy(s, t)
# k = 0.395
# s0 = 0
# s1 = 20
# s_in = 0.5
# t0 = 0
# t1 = 20
# Hxk0 = 30  # кмоль Hxk(t) = Hxk0
# Hxd0 = 50  # кмоль Hxd(t) = Hxd0
# tilde_p = lambda s, t: 0.02
# Fx = lambda t: 86
# a = 738.91
# b = 20
# phi_x = lambda s: a*(s-s_in)**2*np.exp(-b*(s-s_in)) if s < s_in else 0
#
# xfi = (0.222151, 0.14395)
# xi0 = lambda s: 1
# yi0 = lambda s: 0
#
# V = lambda t: 2*np.cos(t/4)             # ЭТО НЕ ТОЧНО !!!!!!!!
# L = lambda t: 3*np.sin(t/2)-30
# L_in = lambda t: Fx(t)*quad(phi_x, s0, s_in)[0]
# dV = lambda t: -1/2*np.sin(t/4)
# dL = lambda t: 3/2*np.cos(t/2)
#
# A11 = lambda s, t: (c1*Fx(t)*phi_x(s)-c1*k*V(t)*tilde_p(s, t)-dL(t))/(L(t)+L_in(t))
# A21 = lambda s, t: c2*k*tilde_p(s, t)
# A12 = lambda s, t: c1*k*V(t)/(L(t)+L_in(t))
# A22 = lambda s, t: - (k*c2-dV(t)/V(t))
# F1 = lambda s, t: c1*Fx(t)*phi_x(s)/(L(t)+L_in(t))
# F2 = lambda s, t: 0
#
# W = lambda t: 2*np.cos(t/4)+3*np.sin(t/2)+56
# G11 = lambda t: (V(t) + W(t))/Hxk0
# G12 = lambda t: -(V(t) + W(t))
# G21 = lambda t: -V(t)/Hxd0
# G22 = lambda t: V(t)/Hxd0
#
#
# rez = sol(s=[s0, s1],
#           c=[c1, c2],
#           fun0=[xi0, yi0],
#           B=[[A11, A12], [A21, A22]],
#           G=[[G11, G12], [G21, G22]],
#           m=10,
#           Time_end=t1)
# print(rez)


