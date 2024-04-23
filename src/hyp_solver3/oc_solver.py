import numpy as np
import scipy

from .problem import HypProblem
from .mesh import Mesh
from .solver import Solver
from typing import List, Callable

ALPHA_HYSTORY = 0
MAX_STEP_ALPHA = 100
MAX_ITERATION = 20


class OCProblem:
    def __init__(self, hyp_problem: HypProblem, mesh: Mesh, xs, ys) -> None:
        self.hyp_problem = hyp_problem
        self.mesh = mesh
        self.solver = Solver()

        self.rezid = []

        self.xs = xs
        self.ys = ys

        nodes = mesh.get_border(type_border="left", sort_t=True)
        t_mesh = [node[0][3] for node in nodes]
        self.t_index = []
        self.t_hash = []
        for i, e in enumerate(t_mesh):
            if not e in self.t_hash:
                self.t_hash.append(e)
                self.t_index.append(i)

    def solve(self):
        self.solver.solve_initial(self.mesh, self.hyp_problem)
        self.solver.solver_center(self.mesh, self.hyp_problem)
        self.solver.solver_final(self.mesh, self.hyp_problem)

    def solve_conj(self):
        self.solver.solver_initial_conj(self.mesh, self.hyp_problem)
        self.solver.solver_center_conj(self.mesh, self.hyp_problem)
        self.solver.solver_final_conj(self.mesh, self.hyp_problem) 

    def J(self, v):
        v_ = [v[ti] for ti in self.t_index]
        self.set_new_control(lambda ti: np.interp(ti, self.t_hash, v_))
        self.solve()

        nodes = self.mesh.get_border(type_border="final", sort_s=True)
        s_final = [node[0][2] for i, node in enumerate(nodes) if i in self.t_index]
        x_final = [node[1][0][0] for i, node in enumerate(nodes)if i in self.t_index]
        y_final = [node[1][0][1] for i, node in enumerate(nodes) if i in self.t_index]

        f = [(xi-self.xs(si))**2+(yi-self.ys(si))**2 for xi, yi, si in zip(x_final, y_final, s_final)]
        return scipy.integrate.trapezoid(f, s_final)


    def set_new_control(self, U):
        G22_0 = lambda t: U(t)
        G21_0 = lambda t: - G22_0(t)
        G11_0 = self.hyp_problem.G11
        G12_0 = self.hyp_problem.G12
        self.hyp_problem.set_G([[G11_0, G12_0], [G21_0, G22_0]])


def solve_oc(hyp_problem: HypProblem, mesh: Mesh, 
             U0: Callable[[float], float], xs: Callable[[float], float], ys: Callable[[float], float], 
             u_lim: List, eps:float=0.000001, debug=True):
    """
    Solve problem:
    J(U) = int_S (x(s, t1)-xs(s))**2 ds + int_S (y(s, t1)-ys(s))**2 ds -> min
    xs, ys - functions

    U0 - start control
    U_lim = [u0, u1] - control limit 
    hyp_problem and mesh are problem
    
    stopping criteria
    int_T H(psi1^k, psi2^k, x^k, y^k, u^k, t)(u^{k+1}-u^k) dt < eps
    """
    oc_problem = OCProblem(hyp_problem, mesh, xs, ys)
    oc_problem.set_new_control(U0)
    for k in range(1, MAX_ITERATION):
        oc_problem.solve()
        oc_problem.solve_conj()
        nodes = mesh.get_border(type_border="left", sort_t=True)
        t = [node[0][3] for node in nodes]
        
        hu = [node[1][1][2]*(node[1][0][0]-node[1][0][1]) for node in nodes]
        
        uk = [hyp_problem.G22(ti) for ti in t]
        new_uk = []
        for uki, hui in zip(uk, hu):
            if hui == 0:
                new_uk.append(uki)
            elif hui > 0:
                new_uk.append(u_lim[1])
            else:
                new_uk.append(u_lim[0])
        

        drez = _get_resid(t, hu, new_uk, uk) 
        oc_problem.rezid.append(drez)
        if debug:
            print("Оптимально" if drez <= eps else "Не оптимально")
            print(f"Невязка составляет: {drez:.6f}")
            print(40*"*")  
        if drez <= eps:
            if debug:
                print(f"Невязка составляет: {drez:.6f}")
                print( f"Решено на {k}-й итерации")
            break
          
        a = _alpha_test(oc_problem, uk, new_uk) 
        uk_a_min = _get_uki(a, uk, new_uk)
        oc_problem.set_new_control(lambda ti: np.interp(ti, t, uk_a_min))
        if debug:
            print("результат минимизации a = ", a)
    return oc_problem


def _get_resid(t_h, h_u_k, uk_new, uk):
    f = [h*(uk_new_i-uk_i) for h, uk_new_i, uk_i in zip(h_u_k, uk_new, uk) ]
    return scipy.integrate.trapezoid(f, t_h) 

def _get_uki(a, uk, uk_new):
    return [ uk_i + a*(uk_new_i - uk_i) for uk_i, uk_new_i in zip(uk, uk_new)]


def _alpha_test(oc_problem:OCProblem, uk, uk_new):
    global ALPHA_HYSTORY
    is_run = True
    is_run_2 = False
    J_old = oc_problem.J(uk)
    a = 2/(2+ALPHA_HYSTORY)
    max_a = ALPHA_HYSTORY + MAX_STEP_ALPHA
    while is_run:
        uk_a = _get_uki(a, uk, uk_new)
        J_new = oc_problem.J(uk_a)
        if is_run_2 and J_old < J_new:
            return a
        if J_old > J_new:
            J_old = J_new
            is_run_2 = True
        ALPHA_HYSTORY += 1
        a =  2/(2+ALPHA_HYSTORY)
        if ALPHA_HYSTORY > max_a:
            return a