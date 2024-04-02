from typing import Tuple
from .problem import HypProblem
import numpy as np
from .mesh import Mesh

class Solver:
    def solve_initial(self, mesh: Mesh, problem):
        for node in mesh.nodes_start_l:
            s = node[2]
            t = node[3]
            mesh.rez_nodes_start_l.append([problem.x0(s), problem.y0(s), 0, 0])
            
        for node in mesh.nodes_start_r:
            s = node[2]
            t = node[3]
            mesh.rez_nodes_start_r.append([problem.x0(s), problem.y0(s), 0, 0])


    def solver_center(self, mesh: Mesh, hyp_problem):
        for node in mesh.nodes_center:
            i = int(node[0])
            j = int(node[1])
            s = node[2]
            t = node[3]      
            # print(f"FACT:{s:.4f}, {t:.4f}", end="")     
            if s==hyp_problem.S0:
                node_l, node_l_rez = mesh.get_s0_node_left(i, j)
                node_r, node_r_rez = mesh.get_center_node_right(i, j) 
            elif s==hyp_problem.S1:
                node_l, node_l_rez = mesh.get_center_node_left(i, j)
                node_r, node_r_rez = mesh.get_s1_node_right(i, j)
            else:
                node_l, node_l_rez = mesh.get_center_node_left(i, j)
                node_r, node_r_rez = mesh.get_center_node_right(i, j)   

            sl, tl = node_l[2], node_l[3]
            sr, tr = node_r[2], node_r[3]
            xl, yl = node_l_rez[0], node_l_rez[1]
            # xl, yl = x_an(sl, tl), y_an(sl, tl)
            xr, yr = node_r_rez[0], node_r_rez[1]
            # xr, yr = x_an(sr, tr), y_an(sr, tr)
            hl = t-tl
            hr = t-tr
            # print(f"{s:.4f}, {t:.4f}| {sl:.4f}, {tl:.4f}, {xl:.4f}, {yl:.4f}, {hl:.4f}| {sr:.4f}, {tr:.4f}, {xr:.4f}, {yr:.4f}, {hr:.4f}")
            # print(f"TRUE:{s:.4f}, {t:.4f}| {x_an(s, t):.4f}, {y_an(s, t):.4f}")
            
            
            if s == hyp_problem.S0:
                x, y = self.left_solve(hyp_problem, s=s, t=t,
                                    sl=sl, tl=tl, xl=xl, yl=yl, hl=hl, 
                                    sr=sr, tr=tr, xr=xr, yr=yr, hr=hr)
            elif s == hyp_problem.S1:
                x, y = self.right_solve(hyp_problem, s=s, t=t,
                                    sl=sl, tl=tl, xl=xl, yl=yl, hl=hl, 
                                    sr=sr, tr=tr, xr=xr, yr=yr, hr=hr)
            else:
                x, y = self.center_solve(hyp_problem, s=s, t=t,
                                    sl=sl, tl=tl, xl=xl, yl=yl, hl=hl, 
                                    sr=sr, tr=tr, xr=xr, yr=yr, hr=hr)
                
            # print(f"| {x:.4f}, {y:.4f}")
            
            mesh.rez_nodes_center.append([x, y])       
             

    def solver_final(self, mesh: Mesh, hyp_problem):
        final_r_nodes = []
        final_l_nodes = []
        for node in mesh.nodes_final_r:
            i, j = node[0], node[1]
            s, t = node[2], node[3]
            
            if s == hyp_problem.S1:
                sr, tr, xr, yr = mesh.get_center_node_stxy(i, j)
                if mesh.is_from_center(i-1, j):
                    s2, t2, x2, y2 = mesh.get_center_node_stxy(i-1, j)
                    sl, tl, xl, yl = mesh.get_stxy_c_3node(sr, tr, xr, yr, s2, t2, x2, y2, s, t, hyp_problem.C2)
                    x, y = self.right_solve(hyp_problem, s=s, t=t,
                                sl=sl, tl=tl, xl=xl, yl=yl, hl=t-tl, 
                                sr=sr, tr=tr, xr=xr, yr=yr, hr=t-tr)
                else:
                    final_r_nodes.append([i, j, sr, tr, xr, yr])
                    x = None; y=None
                mesh.rez_nodes_final_r.append([x, y])
                continue
            sl, tl, xl, yl = mesh.get_center_node_stxy(i, j)
            s1, t1, x1, y1 = mesh.get_center_node_stxy(i+1, j)
            s2, t2, x2, y2 =mesh.get_center_node_stxy(i+1, j+1) if mesh.is_from_center(i+1, j+1) else mesh.get_right_node_stxy(i+1, j)
            sr, tr, xr, yr = mesh.get_stxy_c_3node(s1, t1, x1, y1, s2, t2, x2, y2,s, t, -hyp_problem.C1)
            x, y = self.center_solve(hyp_problem, s=s, t=t,
                                sl=sl, tl=tl, xl=xl, yl=yl, hl=t-tl, 
                                sr=sr, tr=tr, xr=xr, yr=yr, hr=t-tr)
            
            mesh.rez_nodes_final_r.append([x, y]) 

        for node in mesh.nodes_final_l:
            i, j = node[0], node[1]
            s, t = node[2], node[3]
            
            if s == hyp_problem.S0:
                sl, tl, xl, yl = mesh.get_center_node_stxy(i, j)
                if mesh.is_from_center(i, j+1):
                    s2, t2, x2, y2 = mesh.get_center_node_stxy(i, j+1)
                    sr, tr, xr, yr = mesh.get_stxy_c_3node(sl, tl, xl, yl, s2, t2, x2, y2, s, t, -hyp_problem.C1)
                    x, y = self.left_solve(hyp_problem, s=s, t=t,
                                sl=sl, tl=tl, xl=xl, yl=yl, hl=t-tl, 
                                sr=sr, tr=tr, xr=xr, yr=yr, hr=t-tr)
                else:
                    final_l_nodes.append([i, j, sl, tl, xl, yl])
                    x = None; y=None
                mesh.rez_nodes_final_r.append([x, y])
                continue
            sr, tr, xr, yr = mesh.get_center_node_stxy(i, j)
            s1, t1, x1, y1 = mesh.get_center_node_stxy(i, j-1)
            s2, t2, x2, y2 =mesh.get_center_node_stxy(i-1, j-1) if mesh.is_from_center(i-1, j-1) else mesh.get_left_node_stxy(i, j-1)
            sl, tl, xl, yl = mesh.get_stxy_c_3node(s1, t1, x1, y1, s2, t2, x2, y2, s, t, hyp_problem.C2)
            x, y = self.center_solve(hyp_problem, s=s, t=t,
                                sl=sl, tl=tl, xl=xl, yl=yl, hl=t-tl, 
                                sr=sr, tr=tr, xr=xr, yr=yr, hr=t-tr)
             
            mesh.rez_nodes_final_l.append([x, y]) 
        


    def center_solve(self, problem: HypProblem, s:float, t:float,
                    sl:float, tl:float, xl:float, yl:float, hl:float, 
                    sr:float, tr:float, xr:float, yr:float, hr:float) -> Tuple[float, float]:

        B11, B12, B21, B22 = problem.B11, problem.B12, problem.B21, problem.B22,
        F1, F2 = problem.F1, problem.F2
        
        A = [[1 - hr/2*B11(s, t), -hr/2*B12(s, t)],
             [-hl/2*B21(s, t), 1-hl/2*B22(s, t)]]
        b = [xr + hr/2*(B11(sr, tr)*xr+B12(sr, tr)*yr + F1(s, t) + F1(sr, tr)),
             yl + hl/2*(B21(sl, tl)*xl+B22(sl, tl)*yl + F2(s, t) + F2(sl, tl))]
        return np.linalg.solve(A, b)
    

    def left_solve(self, problem: HypProblem, s:float, t:float,
                   sl:float, tl:float, xl:float, yl:float, hl:float, 
                   sr:float, tr:float, xr:float, yr:float, hr:float) -> Tuple[float, float]:
        
        B11 = problem.B11
        B12 = problem.B12
        F1 = problem.F1
        G21 = problem.G21
        G22 = problem.G22

        A = [[1 - hr / 2 * B11(s, t), -hr / 2 * B12(s, t)],
                [-hl / 2 * G21(t), 1 - hl / 2 * G22(t)]]
        b = [xr + hr / 2 * (B11(sr, tr) * xr + B12(sr, tr) * yr + F1(s, t) + F1(sr, tr)),
                yl + hl / 2 * (G21(tl) * xl + G22(tl) * yl)]
        return np.linalg.solve(A, b)
    
    
    def right_solve(self, problem: HypProblem, s:float, t:float,
                   sl:float, tl:float, xl:float, yl:float, hl:float, 
                   sr:float, tr:float, xr:float, yr:float, hr:float) -> Tuple[float, float]:
        
        B21 = problem.B21
        B22 = problem.B22
        F2 = problem.F2
        G12 = problem.G12
        G11 = problem.G11

        A = [[1 - hr / 2 * G11(t), -hr / 2 * G12(t)],
                [-hl / 2 * B21(s, t), 1 - hl / 2 * B22(s, t)]]
        b = [xr + hr / 2 * (G11(tr) * xr + G12(tr) * yr),
                yl + hl / 2 * (B21(sl, tl) * xl + B22(sl, tl) * yl + F2(s, t) + F2(sl, tl))]
        return np.linalg.solve(A, b)