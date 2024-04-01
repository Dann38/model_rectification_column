import unittest
from hyp_solver2.mesh import Mesh
from hyp_solver2.problem import HypProblem
import numpy as np


x0 = lambda s: 0
y0 = lambda s: (3*s+1)/(2*(s+1))

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

PROBLEM_TEST = HypProblem(T=T, S=S, C=C, B=[[B11, B12], [B21, B22]], F=[F1, F2], G=[[G11, G12], [G21, G22]], X0=x0, Y0=y0)
MESH_TEST = Mesh(PROBLEM_TEST, 10)

class NodesTest(unittest.TestCase):
    def test_size_mesh(self):
        mesh = MESH_TEST
        self.assertEqual(len(mesh.nodes_center), 63)
        self.assertEqual(len(mesh.nodes_start_l), 15)
        self.assertEqual(len(mesh.nodes_start_r), 7)

    def test_sorter(self):
        mesh = MESH_TEST
        ar =  [[0, 0, 0, 5], [0, 0, 0, 3], [0, 0, 0, 2], [0, 0, 0, 1]]
        ar_sort = mesh.time_sort(ar).tolist()
        ar.reverse()
       

        for a1, a2 in zip(ar_sort, ar):
            for a1i, a2i in zip(a1, a2):
                self.assertEqual(a1i, a2i)

        ar2_sort = mesh.time_sort(ar).tolist()
        for a1, a2 in zip(ar, ar2_sort):
            for a1i, a2i in zip(a1, a2):
                self.assertEqual(a1i, a2i)

if __name__ == "__main__":
    unittest.main(verbosity=2)