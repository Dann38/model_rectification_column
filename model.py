import numpy as np
from scipy.integrate import quad

# Параметры =====================================
c1 = 36  #м/ч L(s, t)/Hx(s, t)
c2 = 520 #м/ч V(s, t)/Hy(s, t)
k = 0.395
s0 = 0
s1 = 20
s_in = 0.5
t0 = 0
t1 = 20
Hxk0 = 30 #кмоль Hxk(t) = Hxk0
Hxd0 = 50 #кмоль Hxd(t) = Hxd0
tilde_p =  lambda s, t: 0.02
Fx = lambda t: 86
a = 738.91
b = 20
phi_x = lambda s: a*(s-s_in)**2*np.exp(-b*(s-s_in)) if s < s_in else 0

xfi = (0.222151, 0.14395)
xi0 = lambda s: (1, 0)
yi0 = lambda s: (1, 0)

V = lambda t: 2*np.cos(t/4)             # ЭТО НЕ ТОЧНО !!!!!!!!
L = lambda t: 3*np.sin(t/2)-30
L_in = lambda t: Fx(t)*quad(phi_x, s0, s_in)[0]
dV = lambda t: -1/2*np.sim(t/4)
dL = lambda t: 3/2*np.cos(t/2)

A11 = lambda s, t: (c1*Fx(t)*phi_x(s)-c1*k*V(t)*tilde_p(s, t)-dL(t))/(L(t)+L_in(t))
A21 = lambda s, t: c2*k*tilde_p(s, t)
A12 = lambda s, t: c1*k*V(t)/(L(t)+L_in(t))
A22 = lambda s, t: - (k*c2-dV(t)/V(t))
F1 = lambda s, t: c1*Fx(t)*phi_x(s)/(L(t)+L_in(t))
F2 = lambda s, t: 0


G11 = lambda t, W: (V(t) + W(t))/Hxk0
G12 = lambda t, W: -(V(t) + W(t))
G21 = lambda t, W: -V(t)/Hxd0
G22 = lambda t, W: V(t)/Hxd0


