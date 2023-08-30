import matplotlib.pyplot as plt
import numpy as np
from abc import ABC, abstractmethod
C1 = 1
C2 = 3

T0 = 0
T1 = 4
S0 = 0
S1 = 2

B11 = lambda s, t: 1
B12 = lambda s, t: - 1/np.exp(s)
B21 = lambda s, t: np.exp(s)
B22 = lambda s, t: C2


G11 = lambda t:  np.sin(t) / (np.exp(S1) * (np.sin(t) + 2*np.cos(t)) - S1**2 - np.cos(t))
G12 = lambda t: - G11(t)
G21 = lambda t: (2*np.sin(t) - np.cos(t)) / (np.sin(t) + np.cos(t))
G22 = lambda t: - G21(t)

G = [[G11, G12],
     [G21, G22]]

F1 = lambda s, t: -2*C1*s-s**2+np.cos(t)
F2 = lambda s, t: -2*np.exp(s)*np.sin(t)-np.exp(s)*s**2
F = [F1, F2]

x0 = lambda s: s**2+1
y0 = lambda s: 2*np.exp(s)


COUNT_NODE = 60

x_an = lambda s, t: s**2+np.cos(t)
y_an = lambda s, t: np.exp(s)*(np.sin(t)+2*np.cos(t))


class Node(ABC):
    def __init__(self, s, t, left=None, right=None):
        self.s = s
        self.t = t
        self.x = None
        self.y = None
        self.left = left
        self.right = right
        self.is_resolved = False

    def __str__(self):

        # return f"(s: {self.s:.2f}, t: {self.t:.2f})"
        return f"(s: {self.s:.2f}, t: {self.t:.2f}) => x: {self.x: .2f}, y: {self.y: .2f} " \
               f"\t ан.реш x:{x_an(self.s, self.t): .2f}, y:{y_an(self.s, self.t):.2f}"
               # f"({self.left.s:.2f}, {self.left.t:.2f}), ({self.right.s:.2f}, {self.right.t:.2f})"\
               # f"{type(self)}"
               # f"Node left: ({self.left.s: .2f}, {self.left.t: .2f})]"
               # f" [ Node right: ({self.right.s: .2f},{self.right.s: .2f})" \

    @abstractmethod
    def solver(self):
        pass

    def get_inf(self):
        return [self.s, self.t, self.x, self.y]

    def get_old_point(self):

        if not self.left.is_resolved:
            self.left.solver()

        sl, tl, xl, yl = self.left.get_inf()
        hl = self.t - tl

        if not self.right.is_resolved:
            self.right.solver()

        sr, tr, xr, yr = self.right.get_inf()
        hr = self.t - tr
        return (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr)

    def solver_node(self):
        (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = self.get_old_point()

        A = [[1 - hr/2*B11(self.s, self.t), -hr/2*B12(self.s, self.t)],
             [-hl/2*B21(self.s, self.t), 1-hl/2*B22(self.s, self.t)]]
        b = [xr + hr/2*(B11(sr, tr)*xr+B12(sr, tr)*yr + F1(self.s, self.t) + F1(sr, tr)),
             yl + hl/2*(B21(sl, tl)*xl+B22(sl, tl)*yl + F2(self.s, self.t) + F2(sl, tl))]
        self.x, self.y = np.linalg.solve(A, b)


class NodeStart(Node):
    def create_parents_node(self, c, t0):
        if self.left is None:
            self.left = NodeStart(c[1]*(t0-self.t)+self.s, t0)
            self.left.solver()

        if self.right is None:
            s = c[0]*(self.t-t0)+self.s

            self.right = NodeStart(s, t0)
            if s > S1:
                s1t0_node = NodeStart(S1, T0)
                s1t0_node.solver()

                temp_node = NodeRight(self.s, self.t, self.left, s1t0_node)
                temp_node.solver()
                self.x = temp_node.x
                self.y = temp_node.y
                self.is_resolved = True
            else:
                self.right.solver()

    def solver(self):

        if self.t == T0:
            self.x = x0(self.s)
            self.y = y0(self.s)
        elif self.right is None or self.left is None:
            self.create_parents_node([C1, C2], T0)
            if not self.is_resolved:
                self.solver_node()
        else:
            self.solver_node()
        self.is_resolved = True


class NodeLeft(Node):
    def solver(self):
        if self.t == T0:
            self.x = x0(self.s)
            self.y = y0(self.s)
        else:
            (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = self.get_old_point()

            A = [[1 - hr / 2 * B11(self.s, self.t), -hr / 2 * B12(self.s, self.t)],
                 [-hl / 2 * G21(self.t), 1 - hl / 2 * G22(self.t)]]
            b = [xr + hr / 2 * (B11(sr, tr) * xr + B12(sr, tr) * yr + F1(self.s, self.t) + F1(sr, tr)),
                 yl + hl / 2 * (G21(tl) * xl + G22(tl) * yl)]
            self.x, self.y = np.linalg.solve(A, b)
            # self.x = x_an(self.s, self.t)
            # self.y = y_an(self.s, self.t)
        self.is_resolved = True


class NodeRight(Node):
    def solver(self):
        if self.x is not None and self.y is not None:
            pass
        else:
            (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = self.get_old_point()

            A = [[1 - hr / 2 * G11(self.t), -hr / 2 * G12(self.t)],
                 [-hl / 2 * B21(self.s, self.t), 1 - hl / 2 * B22(self.s, self.t)]]
            b = [xr + hr / 2 * (G11(tr) * xr + G12(tr) * yr),
                 yl + hl / 2 * (B21(sl, tl) * xl + B22(sl, tl) * yl + F2(self.s, self.t) + F2(sl, tl))]
            self.x, self.y = np.linalg.solve(A, b)
            # self.x = x_an(self.s, self.t)
            # self.y = y_an(self.s, self.t)
        self.is_resolved = True


class NodeFinish(Node):
    def solver(self):
        self.solver_node()
        self.is_resolved = True


class NoneCenter(Node):
    def solver(self):
        self.solver_node()
        self.is_resolved = True

