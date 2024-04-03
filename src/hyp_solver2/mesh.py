from .problem import HypProblem
import numpy as np
from typing import Tuple

class Mesh:
    def __init__(self, problem: HypProblem, m: int) -> None:
        m_eps = 6
        eps = 0.1**m_eps
        M = m*2
        self.C1 = problem.C1
        self.C2 = problem.C2
        S0 = problem.S0
        S1 = problem.S1
        T0 = problem.T0
        T1 = problem.T1

        Sc = (S1+S0)/2
        Tc = (T1+T0)/2

        DeltaS = (S1-S0)/M
        Delta1T = 1/self.C1*DeltaS 
        Delta2T = 1/self.C2*DeltaS
        DeltaT = Delta1T+Delta2T
        self.dS = DeltaS
        self.d1T = Delta1T
        self.d2T = Delta2T
        
        self.__T_min_ij = (2*T0-(T1-T0))/(2*DeltaS)
        self.__T_max_ij = (2*T1-(T1-T0))/(2*DeltaS)
        self.__S_min_ij = (S0-Sc)/DeltaS
        self.__S_max_ij = (S1-Sc)/DeltaS


        min_j = int((self.__S_max_ij*self.C1+self.__T_max_ij)/(1/self.C1+1/self.C2))
        min_i = int(self.__S_max_ij+min_j) 



        self.__L = np.array([[DeltaS,   DeltaS], 
                    [-Delta1T, Delta2T]])
        self.__v = np.array([Sc, Tc])

        Ind = np.array([[[i, j] for j in range(-min_j, min_j+1)] for i in range(-min_i, min_i+1)])

        Ind_cord = np.array([[(self.__L@vec+self.__v).round(m_eps) for vec in row] for row in Ind])

        nodes_center = []
        nodes_start_l = []
        nodes_start_r = []
        nodes_final_l = []
        nodes_final_r = []

        for row_cord, row_ind in zip(Ind_cord, Ind):
            for vec, inds in zip(row_cord, row_ind) :
                if self.is_from_center(inds[0], inds[1]):
                    nodes_center.append([inds[0], inds[1], vec[0], vec[1]])
                    if vec[0] == S0:
                        if not self.is_from_center(inds[0]+1, inds[1]-1): 
                            nodes_start_l.append([inds[0], inds[1], S0, T0])
                        if not self.is_from_center(inds[0]+1, inds[1]):
                            si = vec[0] + (vec[1]-T0)*self.C1
                            nodes_start_r.append([inds[0], inds[1], si, T0])

                        if not self.is_from_center(inds[0], inds[1]+1):
                            si = vec[0] + (T1 - vec[1])*self.C2 
                            nodes_final_r.append([inds[0], inds[1], si, T1])
                        if not self.is_from_center(inds[0]-1, inds[1]+1):
                            nodes_final_l.append([inds[0], inds[1], S0, T1])
                    elif vec[0] == S1:
                        if not self.is_from_center(inds[0], inds[1]-1):
                            si = vec[0] - (vec[1]-T0)*self.C2 
                            nodes_start_l.append([inds[0], inds[1], si, T0])
                        if not self.is_from_center(inds[0]+1, inds[1]-1):
                            nodes_start_r.append([inds[0], inds[1], S1, T0])

                        if not self.is_from_center(inds[0]-1, inds[1]+1):
                            nodes_final_r.append([inds[0], inds[1], S1, T1])
                        if not self.is_from_center(inds[0]-1, inds[1]):
                            si = vec[0] - (T1 - vec[1])*self.C1
                            nodes_final_l.append([inds[0], inds[1], si, T1])
                    else:
                        if not self.is_from_center(inds[0], inds[1]-1):
                            si = vec[0] - (vec[1]-T0)*self.C2 
                            nodes_start_l.append([inds[0], inds[1], si, T0])
                        if not self.is_from_center(inds[0]+1, inds[1]):
                            si = vec[0] + (vec[1]-T0)*self.C1
                            nodes_start_r.append([inds[0], inds[1], si, T0])

                        if not self.is_from_center(inds[0], inds[1]+1):
                            si = vec[0] + (T1 - vec[1])*self.C2 
                            nodes_final_r.append([inds[0], inds[1], si, T1])
                        if not self.is_from_center(inds[0]-1, inds[1]):
                            si = vec[0] - (T1 - vec[1])*self.C1
                            nodes_final_l.append([inds[0], inds[1], si, T1])


        self.nodes_center = self.time_sort(nodes_center)
        self.nodes_center_dict = self.get_dict(self.nodes_center)

        self.nodes_start_l = self.time_sort(nodes_start_l)
        self.nodes_start_l_dict = self.get_dict(self.nodes_start_l)

        self.nodes_start_r = self.time_sort(nodes_start_r)
        self.nodes_start_r_dict = self.get_dict(self.nodes_start_r)

        self.nodes_final_l = self.time_sort(nodes_final_l)
        self.nodes_final_l_dict = self.get_dict(self.nodes_final_l)

        self.nodes_final_r = self.time_sort(nodes_final_r)
        self.nodes_final_r_dict = self.get_dict(self.nodes_final_r)


        self.rez_nodes_start_l = [[[None, None], [None, None, None]] for i in self.nodes_start_l]
        self.rez_nodes_start_r =  [[[None, None], [None, None, None]] for i in self.nodes_start_r]
        self.rez_nodes_center =  [[[None, None], [None, None, None]] for i in self.nodes_center]
        self.rez_nodes_final_l =  [[[None, None], [None, None, None]] for i in self.nodes_final_l]
        self.rez_nodes_final_r =  [[[None, None], [None, None, None]] for i in self.nodes_final_r]

    def time_sort(self, array):
        ar = np.array(array)
        indexs = [(self.__L @ a[:2]+self.__v)[1] for a in ar]
        index = np.argsort(indexs)
        return ar[index, :]

    def get_dict(self, array):
        nodes_dict = dict()
        
        for index, node in enumerate(array):
            i = int(node[0])
            j = int(node[1])
            if i in nodes_dict.keys():
                nodes_dict[i][j] = index
            else:
                nodes_dict[i] = {j: index}
        return nodes_dict
    
    def get_from_center(self, i, j):
        return self.nodes_center[self.nodes_center_dict[i][j]]

    # def is_from_center(self, i, j):
    #     return i in self.nodes_center_dict.keys() and j in self.nodes_center_dict[i].keys()
    
    def is_from_center(self, i, j):
        t_ind = 1/self.C2*j-1/self.C1*i
        s_ind = i+j
        return (self.__T_min_ij <= t_ind) and (self.__T_max_ij >= t_ind) and \
               (self.__S_min_ij <= s_ind) and (self.__S_max_ij >= s_ind)
    
    def get_center_node(self, i, j):
        node_index = self.nodes_center_dict[i][j]
        node = self.nodes_center[node_index]
        node_rez = self.rez_nodes_center[node_index]
        return node, node_rez
    
    def get_center_node_stxy(self, i, j):
        node, node_rez = self.get_center_node(i, j)
        s, t = node[2], node[3]
        x, y = node_rez[0][0], node_rez[0][1]
        return s, t, x, y

    def get_left_node(self, i, j):
        node_l_index = self.nodes_final_l_dict[i][j]
        node_l = self.nodes_final_l[node_l_index]
        node_l_rez = self.rez_nodes_final_l[node_l_index]
        return node_l, node_l_rez

    def get_left_node_stxy(self, i, j):
        node, node_rez = self.get_left_node(i, j)
        s, t = node[2], node[3]
        x, y = node_rez[0][0], node_rez[0][1]
        return s, t, x, y

    def get_right_node(self, i, j):
        node_r_index = self.nodes_final_r_dict[i][j]
        node_r = self.nodes_final_r[node_r_index]
        node_r_rez = self.rez_nodes_final_r[node_r_index]
        return node_r, node_r_rez
    
    def get_right_node_stxy(self, i, j):
        node, node_rez = self.get_right_node(i, j)
        s, t = node[2], node[3]
        x, y = node_rez[0][0], node_rez[0][1]
        return s, t, x, y

    def get_center_node_left(self, i, j):
        if self.is_from_center(i ,j-1):
            node_l_index = self.nodes_center_dict[i][j-1]
            node_l = self.nodes_center[node_l_index]
            node_l_rez = self.rez_nodes_center[node_l_index]
        else:
            node_l_index = self.nodes_start_l_dict[i][j]
            node_l = self.nodes_start_l[node_l_index]
            node_l_rez = self.rez_nodes_start_l[node_l_index]

        return node_l, node_l_rez
    
    def get_center_node_left_stxy(self, i, j):
        node, node_rez = self.get_center_node_left( i, j)
        s, t = node[2], node[3]
        x, y = node_rez[0][0], node_rez[0][1]
        return s, t, x, y


    def get_center_node_right(self, i, j):
        if self.is_from_center(i+1 ,j):
            node_r_index = self.nodes_center_dict[i+1][j]
            node_r = self.nodes_center[node_r_index]
            node_r_rez = self.rez_nodes_center[node_r_index]
        else:
            node_r_index = self.nodes_start_r_dict[i][j]
            node_r = self.nodes_start_r[node_r_index]
            node_r_rez = self.rez_nodes_start_r[node_r_index]

        return node_r, node_r_rez
    
    def get_center_node_right_stxy(self, i, j):
        node, node_rez = self.get_center_node_right( i, j)
        s, t = node[2], node[3]
        x, y = node_rez[0][0], node_rez[0][1]
        return s, t, x, y

    def get_s0_node_left(self, i, j):
        if self.is_from_center(i+1 ,j-1):
            node_l_index = self.nodes_center_dict[i+1][j-1]
            node_l = self.nodes_center[node_l_index]
            node_l_rez = self.rez_nodes_center[node_l_index]
        else:
            node_l_index = self.nodes_start_l_dict[i][j]
            node_l = self.nodes_start_l[node_l_index]
            node_l_rez = self.rez_nodes_start_l[node_l_index]

        return node_l, node_l_rez
    
    def get_s0_node_left_stxy(self, i, j):
        node, node_rez = self.get_s0_node_left( i, j)
        s, t = node[2], node[3]
        x, y = node_rez[0][0], node_rez[0][1]
        return s, t, x, y
    
    def get_s1_node_right(self, i, j):
        if self.is_from_center(i+1 ,j-1):
            node_r_index = self.nodes_center_dict[i+1][j-1]
            node_r = self.nodes_center[node_r_index]
            node_r_rez = self.rez_nodes_center[node_r_index]
        else:
            node_r_index = self.nodes_start_r_dict[i][j]
            node_r = self.nodes_start_r[node_r_index]
            node_r_rez = self.rez_nodes_start_r[node_r_index]

        return node_r, node_r_rez
    
    def get_s1_node_right_stxy(self, i, j):
        node, node_rez = self.get_s1_node_right( i, j)
        s, t = node[2], node[3]
        x, y = node_rez[0][0], node_rez[0][1]
        return s, t, x, y

    def get_stxy_c_3node(self, s1:float, t1:float, x1:float, y1:float, 
                               s2:float, t2:float, x2:float, y2:float, 
                               s3:float, t3:float, c:float) -> Tuple[float, float, float, float]:
        if t1==t2:
            if t3==t2:
                return s1, t1, x1, y1
            else:
                raise ValueError
        if t1 > t2:
            raise ValueError
        
        c12 = (s1-s2)/(t1-t2)  
        t = 1/(c12-c)*(c12*t1-c*t3 + s3 - s1)
        dt = t-t1
        alpha = dt/(t2-t1)
        s = s1 + alpha*(s2-s1)
        x = x1 + alpha*(x2-x1)
        y = y1 + alpha*(y2-y1)

        return s, t, x, y
    

    # PSI =============================================================
    def get_center_conj_node_stxy(self, i, j):
        node, node_rez = self.get_center_node(i, j)
        s, t = node[2], node[3]
        psi1, psi2, p = node_rez[1][0], node_rez[1][1], node_rez[1][2]
        return s, t, psi1, psi2, p

    def get_final_l_conj_node_stxy(self, i, j):
        node_l_index = self.nodes_final_l_dict[i][j]
        node = self.nodes_final_l[node_l_index]
        node_rez = self.rez_nodes_final_l[node_l_index]
        s, t = node[2], node[3]
        psi1, psi2, p1 = node_rez[1][0], node_rez[1][1], node_rez[1][2]
        return s, t, psi1, psi2, p1

    def get_final_r_conj_node_stxy(self, i, j):
        node_r_index = self.nodes_final_r_dict[i][j]
        node = self.nodes_final_r[node_r_index]
        node_rez = self.rez_nodes_final_r[node_r_index]
        s, t = node[2], node[3]
        psi1, psi2, p2 = node_rez[1][0], node_rez[1][1], node_rez[1][2]
        return s, t, psi1, psi2, p2

    def get_center_node_conj_left_stxy(self, i, j):
        if self.is_from_center(i-1 ,j):
            return self.get_center_conj_node_stxy(i-1 ,j)
        else:
            return self.get_final_l_conj_node_stxy(i, j)

    def get_center_node_conj_right_stxy(self, i, j):
        if self.is_from_center(i ,j+1):
            return self.get_center_conj_node_stxy(i ,j+1)
        else:
            return self.get_final_r_conj_node_stxy(i, j)
    
    def get_s0_node_conj_left_stxy(self, i, j):
        if self.is_from_center(i-1 ,j+1):
            return self.get_center_conj_node_stxy(i-1 ,j+1)
        else:
            return self.get_final_l_conj_node_stxy(i, j)

    def get_s1_node_conj_right_stxy(self, i, j):
        if self.is_from_center(i-1 ,j+1):
            return self.get_center_conj_node_stxy(i-1 ,j+1)
        else:
            return self.get_final_r_conj_node_stxy(i, j)
        
    def get_right_conj_node_stxy(self, i, j):
        node_r_index = self.nodes_start_r_dict[i][j]
        node = self.nodes_start_r[node_r_index]
        node_rez = self.rez_nodes_start_r[node_r_index]
        s, t = node[2], node[3]
        psi1, psi2, p = node_rez[1][0], node_rez[1][1], node_rez[1][2]
        return s, t, psi1, psi2, p

    def get_left_conj_node_stxy(self, i, j):
        node_l_index = self.nodes_start_l_dict[i][j]
        node = self.nodes_start_r[node_l_index]
        node_rez = self.rez_nodes_start_r[node_l_index]
        s, t = node[2], node[3]
        psi1, psi2, p = node_rez[1][0], node_rez[1][1], node_rez[1][2]
        return s, t, psi1, psi2, p