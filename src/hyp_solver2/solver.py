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
            i = node[0]
            j = node[1]
            s = node[2]
            t = node[3]
            
            index_node_l = mesh.nodes_center_dict[i][j]
            node_l = mesh.nodes_center[index_node_l]
            node_l_rez = mesh.rez_nodes_center[index_node_l]
            
            index_node_r1 = mesh.nodes_center_dict[i+1][j] if s!= hyp_problem.S1 else mesh.nodes_center_dict[i+1][j-1]
            node_r1 = mesh.nodes_center[index_node_r1]
            node_r1_rez = mesh.rez_nodes_center[index_node_r1]
            dt =t-node_l[3]
            if mesh.is_from_center(i+1, j+1):
                index_node_r2 = mesh.nodes_center_dict[i+1][j+1]
                node_r2 = mesh.nodes_center[index_node_r2]
                node_r2_rez = mesh.rez_nodes_center[index_node_r2]
            elif s == hyp_problem.S1:                                  # ERROR
                index_node_r2 = mesh.nodes_center_dict[i+2][j-2]
                node_r2 = mesh.nodes_center[index_node_r2]
                node_r2_rez = mesh.rez_nodes_center[index_node_r2]
            else:
                index_node_r2 = mesh.nodes_final_r_dict[i+1][j]
                node_r2 = mesh.nodes_final_r[index_node_r2]
                node_r2_rez = mesh.rez_nodes_final_r[index_node_r2]
            alpha = dt/(node_r2[3]-node_r1[3])

            sl, tl = node_l[2], node_l[3]
            xl, yl = node_l_rez[0], node_l_rez[1]
            hl = t-tl


            sr = node_r1[2]+alpha*(node_r2[2]-node_r1[2])
            if sr >= hyp_problem.S1:
                tr = node_r1[3]+alpha*(node_r2[3]-node_r1[3]) 
                xr = node_r1_rez[0]+alpha*(node_r2_rez[0]-node_r1_rez[0]) 
                yr = node_r1_rez[1]+alpha*(node_r2_rez[1]-node_r1_rez[1])
                sr = hyp_problem.S1
                hr = t-tr
                x, y = self.right_solve(hyp_problem, s=s, t=t,
                                    sl=sl, tl=tl, xl=xl, yl=yl, hl=hl, 
                                    sr=sr, tr=tr, xr=xr, yr=yr, hr=hr)
            else:
                tr = node_r1[3]+alpha*(node_r2[3]-node_r1[3])
                xr = node_r1_rez[0]+alpha*(node_r2_rez[0]-node_r1_rez[0])
                yr = node_r1_rez[1]+alpha*(node_r2_rez[1]-node_r1_rez[1])
                hr = t-tr
                x, y = self.center_solve(hyp_problem, s=s, t=t,
                                    sl=sl, tl=tl, xl=xl, yl=yl, hl=hl, 
                                    sr=sr, tr=tr, xr=xr, yr=yr, hr=hr)
            
            mesh.rez_nodes_final_r.append([x, y]) 

        for node in mesh.nodes_final_l:
            i = node[0]
            j = node[1]
            s = node[2]
            t = node[3]
            
            index_node_r = mesh.nodes_center_dict[i][j]
            node_r = mesh.nodes_center[index_node_r]
            node_r_rez = mesh.rez_nodes_center[index_node_r]
            
            index_node_l1 = mesh.nodes_center_dict[i][j-1] if s!=hyp_problem.S0 else  mesh.nodes_center_dict[i+1][j-1]
            node_l1 = mesh.nodes_center[index_node_l1]
            node_l1_rez = mesh.rez_nodes_center[index_node_l1]
            dt =t-node_r[3]
            if mesh.is_from_center(i-1, j-1):
                index_node_l2 = mesh.nodes_center_dict[i-1][j-1]
                node_l2 = mesh.nodes_center[index_node_l2]
                node_l2_rez = mesh.rez_nodes_center[index_node_l2]
            elif s==hyp_problem.S0:                                  # ERROR
                index_node_l2 = mesh.nodes_center_dict[i+2][j-2] 
                node_l2 = mesh.nodes_center[index_node_l2]
                node_l2_rez = mesh.rez_nodes_center[index_node_l2]
            else:
                index_node_l2 = mesh.nodes_final_l_dict[i][j-1]
                node_l2 = mesh.nodes_final_l[index_node_l2]
                node_l2_rez = mesh.rez_nodes_final_l[index_node_l2]
            alpha = dt/(node_l2[3]-node_l1[3])

            sr, tr = node_r[2], node_r[3]
            xr, yr = node_r_rez[0], node_r_rez[1]
            hr = t-tr

            sl = node_l1[2]+alpha*(node_l2[2]-node_l1[2])
            if sl <= hyp_problem.S0:
                tl = node_l1[3]+alpha*(node_l2[3]-node_l1[3])
                xl = node_l1_rez[0]+alpha*(node_l2_rez[0]-node_l1_rez[0])
                yl = node_l1_rez[1]+alpha*(node_l2_rez[1]-node_l1_rez[1])
                hl = t-tl
                x, y = self.left_solve(hyp_problem, s=s, t=t,
                                    sl=sl, tl=tl, xl=xl, yl=yl, hl=hl, 
                                    sr=sr, tr=tr, xr=xr, yr=yr, hr=hr)
            else:
                tl = node_l1[3]+alpha*(node_l2[3]-node_l1[3])
                xl = node_l1_rez[0]+alpha*(node_l2_rez[0]-node_l1_rez[0])
                yl = node_l1_rez[1]+alpha*(node_l2_rez[1]-node_l1_rez[1])
                hl = t-tl
                x, y = self.center_solve(hyp_problem, s=s, t=t,
                                    sl=sl, tl=tl, xl=xl, yl=yl, hl=hl, 
                                    sr=sr, tr=tr, xr=xr, yr=yr, hr=hr)
           
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