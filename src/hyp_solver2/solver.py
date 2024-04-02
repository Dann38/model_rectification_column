from typing import Tuple
from .problem import HypProblem
import numpy as np
from .mesh import Mesh

class Solver:
    def solve_initial(self, mesh: Mesh, problem):
        for node in mesh.nodes_start_l:
            i, j, s, t = node[0], node[1], node[2], node[3]
            node_index = mesh.nodes_start_l_dict[i][j]
            mesh.rez_nodes_start_l[node_index][0][0]=problem.x0(s) 
            mesh.rez_nodes_start_l[node_index][0][1]=problem.y0(s)
            
        for node in mesh.nodes_start_r:
            i, j, s, t = node[0], node[1], node[2], node[3]
            node_index = mesh.nodes_start_r_dict[i][j]
            mesh.rez_nodes_start_r[node_index][0][0]=problem.x0(s) 
            mesh.rez_nodes_start_r[node_index][0][1]=problem.y0(s)


    def solver_center(self, mesh: Mesh, hyp_problem):
        for node in mesh.nodes_center:
            i, j, s, t = node[0], node[1], node[2], node[3]

            if s==hyp_problem.S0:
                sl, tl, xl, yl = mesh.get_s0_node_left_stxy(i, j)
                sr, tr, xr, yr = mesh.get_center_node_right_stxy(i, j) 
            elif s==hyp_problem.S1:
                sl, tl, xl, yl = mesh.get_center_node_left_stxy(i, j)
                sr, tr, xr, yr = mesh.get_s1_node_right_stxy(i, j)
            else:
                sl, tl, xl, yl = mesh.get_center_node_left_stxy(i, j)
                sr, tr, xr, yr = mesh.get_center_node_right_stxy(i, j)   
            
            if s == hyp_problem.S0:
                x, y = self.left_solve(hyp_problem, s=s, t=t,
                                    sl=sl, tl=tl, xl=xl, yl=yl, hl=t-tl, 
                                    sr=sr, tr=tr, xr=xr, yr=yr, hr=t-tr)
            elif s == hyp_problem.S1:
                x, y = self.right_solve(hyp_problem, s=s, t=t,
                                    sl=sl, tl=tl, xl=xl, yl=yl, hl=t-tl, 
                                    sr=sr, tr=tr, xr=xr, yr=yr, hr=t-tr)
            else:
                x, y = self.center_solve(hyp_problem, s=s, t=t,
                                    sl=sl, tl=tl, xl=xl, yl=yl, hl=t-tl, 
                                    sr=sr, tr=tr, xr=xr, yr=yr, hr=t-tr)
            
            node_index = mesh.nodes_center_dict[i][j]
            mesh.rez_nodes_center[node_index][0][0]=x
            mesh.rez_nodes_center[node_index][0][1]=y      
             

    def solver_final(self, mesh: Mesh, hyp_problem):
        final_r_nodes = []
        for node in mesh.nodes_final_r:
            i, j, s, t = node[0], node[1], node[2], node[3]
            
            if s == hyp_problem.S1:
                sr, tr, xr, yr = mesh.get_center_node_stxy(i, j)
                if mesh.is_from_center(i-1, j):
                    s2, t2, x2, y2 = mesh.get_center_node_stxy(i-1, j)
                    sl, tl, xl, yl = mesh.get_stxy_c_3node(sr, tr, xr, yr, s2, t2, x2, y2, s, t, hyp_problem.C2)
                    x, y = self.right_solve(hyp_problem, s=s, t=t,
                                sl=sl, tl=tl, xl=xl, yl=yl, hl=t-tl, 
                                sr=sr, tr=tr, xr=xr, yr=yr, hr=t-tr)
                    node_index = mesh.nodes_final_r_dict[i][j]
                    mesh.rez_nodes_final_r[node_index][0][0]=x
                    mesh.rez_nodes_final_r[node_index][0][1]=y 
                else:
                    final_r_nodes.append([i, j, sr, tr, xr, yr])
                continue
            sl, tl, xl, yl = mesh.get_center_node_stxy(i, j)
            s1, t1, x1, y1 = mesh.get_center_node_stxy(i+1, j)
            s2, t2, x2, y2 =mesh.get_center_node_stxy(i+1, j+1) if mesh.is_from_center(i+1, j+1) else mesh.get_right_node_stxy(i+1, j)
            sr, tr, xr, yr = mesh.get_stxy_c_3node(s1, t1, x1, y1, s2, t2, x2, y2,s, t, -hyp_problem.C1)
            x, y = self.center_solve(hyp_problem, s=s, t=t,
                                sl=sl, tl=tl, xl=xl, yl=yl, hl=t-tl, 
                                sr=sr, tr=tr, xr=xr, yr=yr, hr=t-tr)
            
            node_index = mesh.nodes_final_r_dict[i][j]
            mesh.rez_nodes_final_r[node_index][0][0]=x
            mesh.rez_nodes_final_r[node_index][0][1]=y 

        for node in mesh.nodes_final_l:
            i, j, s, t = node[0], node[1], node[2], node[3]
            
            if s == hyp_problem.S0:
                sl, tl, xl, yl = mesh.get_center_node_stxy(i, j)
                s2, t2, x2, y2 = mesh.get_center_node_stxy(i, j+1) if mesh.is_from_center(i, j+1) else mesh.get_right_node_stxy(i, j)
                sr, tr, xr, yr = mesh.get_stxy_c_3node(sl, tl, xl, yl, s2, t2, x2, y2, s, t, -hyp_problem.C1)
                x, y = self.left_solve(hyp_problem, s=s, t=t,
                                sl=sl, tl=tl, xl=xl, yl=yl, hl=t-tl, 
                                sr=sr, tr=tr, xr=xr, yr=yr, hr=t-tr)
                node_index = mesh.nodes_final_l_dict[i][j]
                mesh.rez_nodes_final_l[node_index][0][0]=x
                mesh.rez_nodes_final_l[node_index][0][1]=y 
                continue
            sr, tr, xr, yr = mesh.get_center_node_stxy(i, j)
            s1, t1, x1, y1 = mesh.get_center_node_stxy(i, j-1)
            s2, t2, x2, y2 =mesh.get_center_node_stxy(i-1, j-1) if mesh.is_from_center(i-1, j-1) else mesh.get_left_node_stxy(i, j-1)
            sl, tl, xl, yl = mesh.get_stxy_c_3node(s1, t1, x1, y1, s2, t2, x2, y2, s, t, hyp_problem.C2)
            
            x, y = self.center_solve(hyp_problem, s=s, t=t,
                                sl=sl, tl=tl, xl=xl, yl=yl, hl=t-tl, 
                                sr=sr, tr=tr, xr=xr, yr=yr, hr=t-tr)
            
            node_index = mesh.nodes_final_l_dict[i][j]
            mesh.rez_nodes_final_l[node_index][0][0]=x
            mesh.rez_nodes_final_l[node_index][0][1]=y  
        
        for i, j, sr, tr, xr, yr in final_r_nodes:
            s, t, _, _ = mesh.get_right_node_stxy(i, j)
            s2, t2, x2, y2 = mesh.get_left_node_stxy(i, j)
            sl, tl, xl, yl = mesh.get_stxy_c_3node(sr, tr, xr, yr, s2, t2, x2, y2, s, t, hyp_problem.C2)
            x, y = self.right_solve(hyp_problem, s=s, t=t,
                        sl=sl, tl=tl, xl=xl, yl=yl, hl=t-tl, 
                                sr=sr, tr=tr, xr=xr, yr=yr, hr=t-tr)
            node_index = mesh.nodes_final_r_dict[i][j]
            mesh.rez_nodes_final_r[node_index][0][0]=x
            mesh.rez_nodes_final_r[node_index][0][1]=y 


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