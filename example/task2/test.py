from  hyp_solver.problem import HypProblem, ConjHypProblem
import numpy as np
from hyp_solver.mesh import Mesh
from hyp_solver.solve_method import SolveMethod
from hyp_solver.solver import HypSolver
from hyp_solver.ploter import Ploter
import matplotlib.pyplot as plt 

x_an = lambda s, t: 2*(s**2+1)*np.sin(t)
y_an = lambda s, t: (3*s+1)/(2*(s+1))*np.cos(t)
T = [0, 0.4]
S = [0, 1]
C = [1, 2]

B11 = lambda s, t: -2*s
B12 = lambda s, t: 4/3*s**2
B21 = lambda s, t: - (3*s+1)/(4*(s**2+1)*(s+1))
B22 = lambda s, t: 4/((s+1)*(3*s+1))

G11 = lambda t: - 4*np.cos(t) / (np.cos(t)-4*np.sin(t))
G12 = lambda t: - G11(t)
G21 = lambda t: (np.sin(t)) / (np.cos(t) - 4* np.sin(t))
G22 = lambda t: - G21(t)

F1 = lambda s, t: 4*s**3*np.sin(t)+2*np.cos(t)+4*s**2*np.cos(t)/(3*(s+1))
F2 = lambda s, t: 0

x0 = lambda s: x_an(S[1]-s, T[1])
y0 = lambda s: y_an(S[1]-s, T[1])

hyp_problem = HypProblem(T=T, S=S, C=C, B=[[B11, B12], [B21, B22]], F=[F1, F2], G=[[G11, G12], [G21, G22]], X0=x0, Y0=y0)
conj_hyp_problem = ConjHypProblem(hyp_problem)

mesh = Mesh(conj_hyp_problem, count_node=20)
method = SolveMethod()
solver = HypSolver(mesh, method)

solver.solve(conj_hyp_problem)
ploter = Ploter()

x_an_rez = lambda s, t: x_an(S[1]-s, T[1]-t)
y_an_rez = lambda s, t: y_an(S[1]-s, T[1]-t)

ploter.plot_final(mesh, x_an_rez, y_an_rez)
