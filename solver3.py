import matplotlib.pyplot as plt
import numpy as np
from abc import ABC, abstractmethod
C1 = 1
C2 = 3

x0 = lambda s: np.exp(s)
y0 = lambda s: 0

B11 = lambda s, t: - C1
B12 = lambda s, t: - np.exp(s)/(s+2)
B21 = lambda s, t: (s+2)/np.exp(s)
B22 = lambda s, t: C2/(s+2)

T0 = 0
T1 = 4

COUNT_NODE = 15

x_an = lambda s, t: np.exp(s)*np.cos(t)
y_an = lambda s, t: (s+2)*np.sin(t)


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
               f"\t ан.реш x:{x_an(self.s, self.t): .2f}, y:{y_an(self.s, self.t):.2f}" \
               f"({self.left.s:.2f}, {self.left.t:.2f}), ({self.right.s:.2f}, {self.right.t:.2f})"
               # f"Node left: ({self.left.s: .2f}, {self.left.t: .2f})]"
               # f" [ Node right: ({self.right.s: .2f},{self.right.s: .2f})" \

    @abstractmethod
    def solver(self):
        pass

    def solver_node(self):
        if not self.left.is_resolved:
            self.left.solver()

        xl = self.left.x
        yl = self.left.y
        sl = self.left.s
        tl = self.left.t
        hl = self.t - tl

        if not self.right.is_resolved:
            self.right.solver()

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
    def create_parents_node(self, c, t0):
        if self.left is None:
            self.left = NodeStart(c[1]*(t0-self.t)+self.s, t0)
            self.left.solver()

        if self.right is None:
            self.right = NodeStart(c[0]*(self.t-t0)+self.s, t0)
            self.right.solver()

    def solver(self):
        if self.t == T0:
            self.x = x0(self.s)
            self.y = y0(self.s)
        elif self.right is None or self.left is None:
            self.create_parents_node([C1, C2], T0)
            self.solver_node()
        else:
            self.solver_node()


class NodeLeft(Node):
    def solver(self):
        self.x = x_an(self.s, self.t)
        self.y = y_an(self.s, self.t)
        self.is_resolved = True


class NodeRight(Node):
    def solver(self):
        self.x = x_an(self.s, self.t)
        self.y = y_an(self.s, self.t)
        self.is_resolved = True


class NodeFinish(Node):
    def solver(self):
        pass


class NoneCenter(Node):
    def solver(self):
        self.solver_node()
        self.is_resolved = True


class Mesh:
    def __init__(self):
        self.c1 = C1
        self.c2 = C2
        self.t0 = T0
        self.t1 = T1
        self.s0 = 0
        self.s1 = 2
        self.m = COUNT_NODE
        self.ds, self.h_down, self.h_up = self.create_mesh_s()
        self.start = MeshStart(self.ds, [self.h_down, self.h_up], self.m, self.s0, self.t0)
        inds = self.start.neighbor_is_left()
        self.left = MeshLeft([self.h_down, self.h_up], self.s0, self.t1)
        self.right = MeshRight([self.h_down, self.h_up], self.s1, self.t1)
        self.center = MeshCenter(self.ds, [self.h_down, self.h_up], self.m, inds)
        self.finish = MeshFinish()

        self.result = None

    def create_mesh_s(self):
        ds = (self.s1 - self.s0) / self.m
        h_down = ds / self.c2
        h_up = ds / self.c1
        return ds, h_down, h_up

    def solve(self):
        start_nodes = self.start.get_nodes()

        # print("получение стартовых узлов")
        for node in start_nodes:  # начальные условия
            node.solver()

        # print("получение оставшихся узлов")
        left_nodes = self.left.get_nodes(start_nodes[0])
        right_nodes = self.right.get_nodes(start_nodes[-1])
        center_nodes = self.center.get_nodes(start_nodes[1:-1], left_nodes, right_nodes)


        # print("решение узлов по уровням")
        for i, level_nodes in enumerate(center_nodes):
            left_nodes[i].solver()
            for node in level_nodes:
                node.solver()
            right_nodes[i].solver()
        finish_nodes = self.finish.get_nodes(left_nodes[-1], center_nodes[-1], right_nodes[-1])
        self.result = (start_nodes, left_nodes, right_nodes, center_nodes)

    def plot_mesh(self):
        if self.result is None:
            self.solve()
        start_nodes, left_nodes, right_nodes, center_nodes = self.result

        def plot_node(node, c="g"):
            plt.scatter(node.s, node.t, color=c)
            # if node.left.x is None and node.left.y is None:
            #     plt.scatter(node.s, node.t, color="r")
            # elif node.left.x is None:
            #     plt.scatter(node.s, node.t, color="b")
            # elif node.left.y is None:
            #     plt.scatter(node.s, node.t, color="orange")
            # else:
            #     plt.scatter(node.s, node.t, color="g")

        for node in start_nodes:
            plot_node(node, "r")
        for node in left_nodes:
            plot_node(node)
        for node in right_nodes:
            plot_node(node)
        for level in center_nodes:
            for node in level:
                plot_node(node)
        plt.xlim([self.s0, self.s1])
        plt.ylim([self.t0, self.t1])
        plt.show()

    def get_XY_t1(self):
        if self.result is None:
            self.solve()
        start_nodes, left_nodes, right_nodes, center_nodes = self.result

        def approx(t_u, t_d, f_u, f_d):

            if t_d == t_u:
                return f_u
            return (f_d-f_u)*(T1-t_u)/(t_d-t_u) + f_u

        nodes_up = [left_nodes[-2]] + [c for c in center_nodes[-2]] + [right_nodes[-2]]
        nodes_down = [left_nodes[-3]] + [c for c in center_nodes[-3]] + [right_nodes[-3]]

        x_t1 = [approx(u.t, d.t, u.x, d.x) for u, d in zip(nodes_up, nodes_down)]
        y_t1 = [approx(u.t, d.t, u.y, d.y) for u, d in zip(nodes_up, nodes_down)]
        s_ = [n.s for n in nodes_up]

        return s_, x_t1, y_t1


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
            t_ = (t_ + self.h[0]) % (self.h[0]+self.h[1])
            s_ += self.ds

        inds = self.neighbor_is_left()
        for i, ind in enumerate(inds):
            if ind == 1.:
                nodes[i+1].left = nodes[i]
            else:
                nodes[i].right = nodes[i+1]

        return nodes

    def neighbor_is_left(self):
        ind = np.zeros((self.m))
        t_hist = 0
        t_ = (t_hist + self.h[0]) % (self.h[0] + self.h[1])
        for i in range(self.m):
            if t_ > t_hist:
                ind[i] = 1.
            t_hist = t_
            t_ = (t_hist + self.h[0]) % (self.h[0] + self.h[1])
        return ind


class MeshCenter:
    def __init__(self, ds, h, m, inds):
        self.ds = ds
        self.dt = h[0]+h[1]
        self.h = h
        self.m = m
        self.inds = inds

    def get_nodes(self, start_nodes, left_nodes, right_nodes):
        count_level = len(left_nodes)
        nodes = [[NoneCenter(node.s, node.t+self.dt*i) for node in start_nodes] for i in range(1, count_level)]
        # Стыковка левых точек ----------------------------------------------------------------------------------------
        l_start = 1
        for i, left_node in enumerate(left_nodes[l_start:]):
            nodes[i][0].left = left_node

        for i in range(2, count_level):
            left_nodes[i].right = nodes[i-2][0]
        left_nodes[1].right = start_nodes[0]
        # Стыковка правых точек ---------------------------------------------------------------------------------------

        if self.inds[-1] == 0:
            for i, right_node in enumerate(right_nodes[1:]):
                nodes[i][-1].right = right_node
            for i in range(1, len(right_nodes)-1):
                right_nodes[1+i].left = nodes[i-1][-1]
            right_nodes[1].left = start_nodes[-1]
        else:
            right_nodes[0].left = start_nodes[-1]
            for i in range(1, count_level):
                right_nodes[i].left = nodes[i-1][-1]
            for i in range(count_level-1):
                nodes[i][-1].right = right_nodes[i]



        # Стыковка начальных точек ------------------------------------------------------------------------------------
        for i, node in enumerate(start_nodes[1:]):
            if self.inds[i+1] == 1:
                nodes[0][i].right = node

        # Стыковка соседних точек в слое ------------------------------------------------------------------------------
        for level in range(0, count_level-1): # COUNT_LEVEL - вместе со стартовым
            for i, ind in enumerate(self.inds[1:-1]):
                if ind == 1.:
                    nodes[level][i+1].left = nodes[level][i]
                    if level == 0:
                        nodes[level][i].right = start_nodes[i+1]
                    else:
                        nodes[level][i].right = nodes[level-1][i+1]
                else:
                    nodes[level][i].right = nodes[level][i+1]
                    if level == 0:
                        nodes[level][i + 1].left = start_nodes[i]
                    else:
                        nodes[level][i + 1].left = nodes[level-1][i]

        return nodes


class MeshLeft:
    def __init__(self, h, s0, t1):
        self.h = h
        self.s0 = s0
        self.t1 = t1

    def get_nodes(self, start_node):
        dt = self.h[0]+self.h[1]
        t = np.arange(start_node.t+dt, self.t1, dt)

        start_left_node = NodeLeft(start_node.s, start_node.t, right=start_node.right)
        start_left_node.x = start_node.x
        start_left_node.y = start_node.y

        nodes = [start_left_node]

        nodes = nodes + [NodeLeft(self.s0, ti) for ti in t]
        for i, node in enumerate(nodes[1:]):
            node.left = nodes[i]
        return nodes


class MeshRight:
    def __init__(self, h, s1, t1):
        self.h = h
        self.s1 = s1
        self.t1 = t1


    def get_nodes(self, start_node):
        dt = self.h[0] + self.h[1]
        t = np.arange(start_node.t+dt, self.t1, dt)

        start_right_node = NodeRight(start_node.s, start_node.t, left=start_node.left)
        start_right_node.x = start_node.x
        start_right_node.y = start_node.y

        nodes = [start_right_node]
        nodes = nodes + [NodeRight(self.s1, ti) for ti in t]

        for i, node in enumerate(nodes):
            node.right = nodes[i]
        return nodes


class MeshFinish:
    def get_nodes(self, left_node, center_nodes, right_node):
        print(left_node)
        for node in center_nodes:
            print(node)
        print(right_node)


if __name__ == '__main__':
    mesh = Mesh()
    mesh.solve()
    # mesh.plot_mesh()
    # s, x, y = mesh.get_XY_t1()
    # x_an_ = [x_an(si, T1) for si in s]
    # y_an_ = [y_an(si, T1) for si in s]
    #
    # fig, axs = plt.subplots(nrows=2, ncols=1)
    # axs[0].set_title(f"Численное решение (узлов по t: 232, узлов по s: {COUNT_NODE})")
    # axs[0].plot(s, x_an_)
    # axs[0].plot(s, x, "o")
    # axs[0].grid()
    # axs[0].legend(["аналитическое решение", "численное решение"])
    # axs[0].set_ylabel("$x(s, t_1)$")
    #
    # axs[1].plot(s, y_an_)
    # axs[1].plot(s, y, "o")
    # axs[1].grid()
    # axs[1].legend(["аналитическое решение", "численное решение"])
    # axs[1].set_ylabel("$y(s, t_1)$")
    # axs[1].set_xlabel("$s$")
    #
    # print(np.max(abs(np.array(x_an_)-np.array(x))))
    # print(np.max(abs(np.array(y_an_)-np.array(y))))
    # plt.show()
