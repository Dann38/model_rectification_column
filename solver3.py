import matplotlib.pyplot as plt
import numpy as np

C1 = 1
C2 = 3

x0 = lambda s: np.exp(s)
y0 = lambda s: 0

B11 = lambda s, t: - C1
B12 = lambda s, t: - np.exp(s)/(s+2)
B21 = lambda s, t: (s+2)/np.exp(s)
B22 = lambda s, t: C2/(s+2)

T0 = 0

x_an = lambda s, t: np.exp(s)*np.cos(t)
y_an = lambda s, t: (s+2)*np.sin(t)

class Node:
    def __init__(self, s, t, left=None, right=None):
        self.s = s
        self.t = t
        self.x = None
        self.y = None
        self.left = left
        self.right = right

    def __str__(self):
        return f"(s: {self.s:.2f}, t: {self.t:.2f}) => x: {self.x: .2f}, y: {self.y: .2f} " \
               f"(ан: ({x_an(self.s, self.t): .2f}, {y_an(self.s, self.t):.2f}))"

               # f" [ Node right: ({self.right.s: .2f},{self.right.s: .2f})" \
               # f"Node left: ({self.left.s: .2f},{self.left.s: .2f})]"

    def solver_node(self):

        if self.left.x is None or self.left.y is None:
            print("in left")
            self.left.solver_node()

        xl = self.left.x
        yl = self.left.y
        sl = self.left.s
        tl = self.left.t
        hl = self.t - tl

        if self.right.x is None or self.right.y is None:
            print("in right")
            self.right.solver_node()

        xr = self.right.x
        yr = self.right.y
        sr = self.right.s
        tr = self.right.t
        hr = self.t - tr

        A = [[1 - hr/2*B11(self.s, self.t), -hr/2*B12(self.s, self.t)],
             [-hl/2*B21(self.s, self.t), 1-hl/2*B22(self.s, self.t)]]
        b = [xr + hr/2*(B11(sr, tr)*xr+B12(sr, tr)*yr),
             yl + hl/2*(B21(sl, tl)*xl+B22(sl, tl)*yl)]
        self.x, self.y = np.linalg.solve(A, b)


class NodeStart(Node):
    def get_parents_node(self, c, t0, h):
        if self.t - h[0] <= 0:
            left = NodeStart(c[1]*(t0-self.t)+self.s, t0)
            self.left = left
        else:
            left = self.left

        if self.t - h[1] <= 0:
            right = NodeStart(c[0]*(self.t-t0)+self.s, t0)
            self.right = right
        else:
            right = self.right
        return left, right

    def solver_t0(self, t0):
        if self.t == t0:
            self.x = x0(self.s)
            self.y = y0(self.s)


class Mesh:
    def __init__(self):
        self.c1 = C1
        self.c2 = C2
        self.t0 = T0
        self.t1 = 4
        self.s0 = 0
        self.s1 = 2
        self.m = 30
        self.ds, self.h_down, self.h_up = self.create_mesh_s()
        self.start = MeshStart(self.ds, [self.h_down, self.h_up], self.m, self.s0, self.t0)
        self.center = MeshCenter()
        self.left = MeshLeft()
        self.right = MeshRight()

    def create_mesh_s(self):
        ds = (self.s1 - self.s0) / self.m
        h_down = ds / self.c2
        h_up = ds / self.c1
        return ds, h_down, h_up

    def plot(self):
        nodes = self.start.get_nodes()
        for node in nodes:
            l, r = node.get_parents_node([self.c1, self.c2], self.t0, [self.h_down, self.h_up])
            l.solver_t0(self.t0)
            r.solver_t0(self.t0)
            node.solver_t0(self.t0)

        for node in nodes:
            if node.t != self.t0:
               node.solver_node()
            print(node)

        for node in nodes:
            plt.scatter(node.s, node.t, color="g")
            l, r = node.get_parents_node([self.c1, self.c2], self.t0, [self.h_down, self.h_up])
            plt.scatter(l.s, l.t, color="r", marker=".")
            plt.scatter(r.s, r.t, color="r", marker=".")

        plt.xlim([self.s0, self.s1])
        plt.ylim([self.t0, self.t1])
        plt.show()


class MeshStart:
    def __init__(self, ds, h,  m, s0, t0):
        self.ds = ds
        self.h = h
        self.m = m
        self.s0 = s0
        self.t0 = t0

    def get_nodes(self):

        nodes = []
        t_ = self.t0
        s_ = self.s0
        for i in range(self.m+1):
            nodes.append(NodeStart(s_, t_))
            if i != 0:
                if nodes[-1].t > nodes[-2].t:
                    nodes[-1].left = nodes[-2]
                if nodes[-2].t > nodes[-1].t:
                    nodes[-1].right = nodes[-2]
            t_ = (t_ + self.h[0]) % (self.h[0]+self.h[1])
            s_ += self.ds

        return nodes


class MeshCenter:
    pass


class MeshLeft:
    pass


class MeshRight:
    pass


if __name__ == '__main__':
    mesh = Mesh()
    mesh.plot()
