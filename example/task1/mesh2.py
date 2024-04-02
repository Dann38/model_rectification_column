from hyp_solver2 import HypProblem, Mesh, Solver
import matplotlib.pyplot as plt
import numpy as np

x0 = lambda s: 0
y0 = lambda s: (3*s+1)/(2*(s+1))

x_an = lambda s, t: 2*(s**2+1)*np.sin(t)
y_an = lambda s, t: (3*s+1)/(2*(s+1))*np.cos(t)
T = [0, 0.4]
S = [0, 1]
C = [1,2]
# C = [2.3,0.9]

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

hyp_problem = HypProblem(T=T, S=S, C=C, B=[[B11, B12], [B21, B22]], F=[F1, F2], G=[[G11, G12], [G21, G22]], X0=x0, Y0=y0)

mesh = Mesh(hyp_problem, 10)

for node in mesh.nodes_center:
    i = int(node[0])
    j = int(node[1])
    s = node[2]
    t = node[3]
    if mesh.is_from_center(i, j-1):
        s_, t_, _, _ = mesh.get_center_node_stxy(i, j-1)
        plt.plot([s,s_], [t, t_], "g" )
    if mesh.is_from_center(i+1, j):
        s_, t_, _, _ = mesh.get_center_node_stxy(i+1, j)
        plt.plot([s,s_], [t, t_], "g")
    plt.scatter(x=s, y=t, color="b")

for node in mesh.nodes_final_l:
    i = int(node[0])
    j = int(node[1])
    s = node[2]
    t = node[3]
    
    s_, t_, _, _ = mesh.get_center_node_stxy(i, j)
    plt.plot([s,s_], [t, t_], "g")
    plt.scatter(s, t, color="r", marker=".")

for node in mesh.nodes_final_r:
    i = int(node[0])
    j = int(node[1])
    s = node[2]
    t = node[3]
    
    s_, t_, _, _ = mesh.get_center_node_stxy(i, j)
    plt.plot([s,s_], [t, t_], "g")
    plt.scatter(node[2], node[3], color="r", marker=".")
    
for node in mesh.nodes_start_l:
    i = int(node[0])
    j = int(node[1])
    s = node[2]
    t = node[3]
    node, _ = mesh.get_center_node(i, j)
    plt.plot([s,node[2]], [t, node[3]], "g")
    plt.scatter(s, t, color="r", marker=".")

for node in mesh.nodes_start_r:
    i = int(node[0])
    j = int(node[1])
    s = node[2]
    t = node[3]
    node, _ = mesh.get_center_node(i, j)
    plt.plot([s,node[2]], [t, node[3]], "g")
    plt.scatter(s, t, color="r", marker=".")


plt.plot([hyp_problem.S0, hyp_problem.S1], [hyp_problem.T0, hyp_problem.T0], "r")
plt.plot([hyp_problem.S0, hyp_problem.S1], [hyp_problem.T1, hyp_problem.T1], "r")
plt.plot([hyp_problem.S0, hyp_problem.S0], [hyp_problem.T0, hyp_problem.T1], "r")
plt.plot([hyp_problem.S1, hyp_problem.S1], [hyp_problem.T0, hyp_problem.T1], "r")


plt.xlabel("s")
plt.ylabel("t")
plt.savefig("mesh",format="eps")