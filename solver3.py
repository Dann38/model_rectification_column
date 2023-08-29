import matplotlib.pyplot as plt
import numpy as np
from abc import ABC, abstractmethod
C1 = 1
C2 = 3

T0 = 0
T1 = 4
S0 = 0
S1 = 2

x0 = lambda s: np.exp(s)
y0 = lambda s: 0

B11 = lambda s, t: - C1
B12 = lambda s, t: - np.exp(s)/(s+2)
B21 = lambda s, t: (s+2)/np.exp(s)
B22 = lambda s, t: C2/(s+2)

G11 = lambda t: np.exp(S1)*np.sin(t)/((S1+2)*np.sin(t)-np.exp(S1)*np.cos(t))
G12 = lambda t: - G11(t)
G21 = lambda t: 2*np.cos(t)/(np.cos(t)-2*np.sin(t))
G22 = lambda t: - G21(t)

COUNT_NODE = 24

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
            b = [xr + hr / 2 * (B11(sr, tr) * xr + B12(sr, tr) * yr),
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
                 [-hl / 2 * B21(self.t, self.t), 1 - hl / 2 * B22(self.t, self.t)]]
            b = [xr + hr / 2 * (G11(tr) * xr + G12(tr) * yr),
                 yl + hl / 2 * (B21(tl, tr) * xl + B22(tl, tr) * yl)]
            self.x, self.y = np.linalg.solve(A, b)
            # self.x = x_an(self.s, self.t)
            # self.y = y_an(self.s, self.t)
        self.is_resolved = True


class NodeFinish(Node):
    def solver(self):
        self.solver_node()
        self.is_resolved = True

    # def get_xy_in_t1(self):
    #     self.

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
        self.s0 = S0
        self.s1 = S1
        self.m = COUNT_NODE
        self.ds, self.h_down, self.h_up = self.create_mesh_s()
        self.start = MeshStart(self.ds, [self.h_down, self.h_up], self.m, self.s0, self.t0)
        inds = self.start.neighbor_is_left()
        self.left = MeshLeft([self.h_down, self.h_up], self.s0, self.t1)
        self.right = MeshRight([self.h_down, self.h_up], self.s1, self.t1)
        self.center = MeshCenter(self.ds, [self.h_down, self.h_up], self.m, inds)
        self.finish = MeshFinish([self.h_down, self.h_up], inds)

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
        right_nodes = self.right.get_nodes(start_nodes[-1], len(left_nodes))
        center_nodes = self.center.get_nodes(start_nodes[1:-1], left_nodes, right_nodes)


        # print("решение узлов по уровням")
        for i, level_nodes in enumerate(center_nodes):
            left_nodes[i+1].solver()
            for node in level_nodes:
                node.solver()
            right_nodes[i+1].solver()

        finish_nodes = self.finish.get_nodes(left_nodes[-1], center_nodes[-1], right_nodes[-1], T1)

        for node in finish_nodes:
            node.solver()


        finish_nodes2 = self.finish.get_nodes(finish_nodes[0], finish_nodes[1:-1], finish_nodes[-1],
                                              T1+self.h_down+self.h_up)

        for node in finish_nodes2:
            node.solver()
        self.result = (start_nodes, left_nodes, right_nodes, center_nodes, finish_nodes, finish_nodes2)

    def plot_mesh(self):
        # if self.result is None:
        #     self.solve()
        start_nodes, left_nodes, right_nodes, center_nodes, finish_nodes, finish_nodes2 = self.result

        fig, ax = plt.subplots()


        def onclick(event):
            print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
                  ('double' if event.dblclick else 'single', event.button,
                   event.x, event.y, event.xdata, event.ydata))

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        def plot_node(node, c="g"):
            # plt.scatter(node.s, node.t, color=c)
            s = f"{node.t:.2f}"
            if node.left is None and node.right is None:
                plt.scatter(node.s, node.t, color="r")
                s += f" "
            elif node.left is None:
                plt.scatter(node.s, node.t, color="b")
                s += f" ( , {node.right.t:.2f})"
            elif node.right is None:
                plt.scatter(node.s, node.t, color="orange")
                s += f" ({node.left.t:.2f}, )"
            else:
                plt.scatter(node.s, node.t, color="g")
                s += f" ({node.left.t:.2f}, {node.right.t:.2f})"
            plt.text(node.s, node.t, s, fontsize="xx-small")

        for node in start_nodes:
            plot_node(node)
        for node in left_nodes:
            plot_node(node)
        for node in right_nodes:
            plot_node(node)
        for level in center_nodes:
            for node in level:
                plot_node(node)
        for node in finish_nodes:
            plot_node(node)
        for node in finish_nodes2:
            plot_node(node)

        # x = [node.x for node in left_nodes]
        # y = [node.y for node in left_nodes]
        # t = [node.t for node in left_nodes]
        # plt.plot(t, x)
        # plt.plot(np.array(t), x_an(S0, np.array(t)))

        plt.xlim([self.s0, self.s1])
        plt.ylim([self.t0, self.t1+0.5])
        plt.show()

    def get_XY_t1(self):
        if self.result is None:
            self.solve()
        start_nodes, left_nodes, right_nodes, center_nodes, f1, f2 = self.result

        def approx(t_u, t_d, f_u, f_d):

            if t_d == t_u:
                return f_u
            return (f_d-f_u)*(T1-t_u)/(t_d-t_u) + f_u


        x_t1 = [approx(u.t, d.t, u.x, d.x) for u, d in zip(f1, f2)]
        y_t1 = [approx(u.t, d.t, u.y, d.y) for u, d in zip(f1, f2)]
        s_ = [n.s for n in f1]

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

        # Стыковка начальных точек ------------------------------------------------------------------------------------
        for i, node in enumerate(start_nodes[1:]):
            if self.inds[i + 1] == 1:
                nodes[0][i].right = node

        # Стыковка левых точек ----------------------------------------------------------------------------------------
        for i in range(count_level-1):
            nodes[i][0].left = left_nodes[i+1]

        for i in range(count_level-2):
            # print(f"l: {left_nodes[i+2].s:.2f}, {left_nodes[i+2].t:.2f}   <- {nodes[i][0].s:.2f}, {nodes[i][0].t:.2f}")
            left_nodes[i+2].right = nodes[i][0]

        left_nodes[1].right = start_nodes[0]


        # Стыковка правых точек ---------------------------------------------------------------------------------------

        if self.inds[-1] == 1.:
            print("ok")
            for level in range(0, count_level - 1):  # COUNT_LEVEL - вместе со стартовым
                right_nodes[1+level].left = nodes[level][-1]
                nodes[level][-1].right = right_nodes[level]
            # print(nodes[count_level-1].s, nodes[count_level-1].t)
        else:
            for level in range(0, count_level - 2):
                right_nodes[2+level].left = nodes[level][-1]
                nodes[level][-1].right = right_nodes[level+1]
            nodes[-1][-1].right = right_nodes[-1]
            right_nodes[1].left = start_nodes[-1]

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
        # num = 3
        # print(nodes[-1][num].s, nodes[-1][num].t)
        # print(nodes[-1][num].left.s, nodes[-1][num].left.t, nodes[-1][num].right.s, nodes[-1][num].right.t)
        return nodes


class MeshLeft:
    def __init__(self, h, s0, t1):
        self.h = h
        self.s0 = s0
        self.t1 = t1

    def get_nodes(self, start_node):
        dt = self.h[0]+self.h[1]
        t = np.arange(start_node.t+dt, self.t1-dt, dt)

        start_left_node = NodeLeft(start_node.s, start_node.t)
        start_left_node.is_resolved = True

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


    def get_nodes(self, start_node, count_nodes):
        dt = self.h[0] + self.h[1]
        t = [ start_node.t+i*dt for i in range(1, count_nodes)]

        start_right_node = NodeRight(start_node.s, start_node.t, left=start_node.left)
        start_right_node.is_resolved = True
        start_right_node.x = start_node.x
        start_right_node.y = start_node.y

        nodes = [start_right_node]
        nodes = nodes + [NodeRight(self.s1, ti) for ti in t]

        for i, node in enumerate(nodes[1:]):
            node.right = nodes[i]
        return nodes


class MeshFinish:
    def __init__(self, h, inds):
        self.dt = h[1]+h[0]
        self.inds = inds

    def get_nodes(self, left_node, center_nodes, right_node, t1):
        level_end = [left_node] + center_nodes + [right_node]
        finish_nodes = [NoneCenter(node.s, node.t+self.dt) for node in level_end]

        fin_nodes = []
        for i, f in enumerate(finish_nodes):
            if f.t < t1:
                if i == 0:
                    fin_nodes.append(NodeLeft(f.s, f.t))
                elif i == len(finish_nodes)-1:
                    fin_nodes.append(NodeRight(f.s, f.t))
                else:
                    fin_nodes.append(NodeFinish(f.s, f.t))
            else:
                fin_nodes.append(level_end[i])

        for i in range(0, len(fin_nodes)-1):
            if fin_nodes[i+1].t < fin_nodes[i].t:
                if fin_nodes[i].right is None:
                    fin_nodes[i].right = fin_nodes[i+1]
                if fin_nodes[i+1].left is None:
                    fin_nodes[i+1].left = level_end[i]
            else:
                if fin_nodes[i + 1].left is None:
                    fin_nodes[i+1].left = fin_nodes[i]
                if fin_nodes[i].right is None:
                    fin_nodes[i].right = level_end[i+1]

        # Правая граница для правого края ---------------------------
        if fin_nodes[-1].t == right_node.t:
            fin_nodes[-1] = right_node
        else:
            fin_nodes[-1].right = right_node
        # Левая граница для левого края -----------------------------
        if fin_nodes[1].t > fin_nodes[0].t:
            fin_nodes[0].right = level_end[1]
        else:
            fin_nodes[0].right = fin_nodes[1]
        fin_nodes[0].left = level_end[0]
        # -----------------------------------------------------------

        return fin_nodes


if __name__ == '__main__':
    mesh = Mesh()
    mesh.solve()
    for node in mesh.result[-3][-2]:
        print(node)
    # mesh.plot_mesh()
    s, x, y = mesh.get_XY_t1()
    x_an_ = [x_an(si, T1) for si in s]
    y_an_ = [y_an(si, T1) for si in s]

    fig, axs = plt.subplots(nrows=2, ncols=1)
    axs[0].set_title(f"Численное решение (узлов по s: {COUNT_NODE})")
    axs[0].plot(s, x_an_)
    axs[0].plot(s, x, "o")
    axs[0].grid()
    axs[0].legend(["аналитическое решение", "численное решение"])
    axs[0].set_ylabel("$x(s, t_1)$")

    axs[1].plot(s, y_an_)
    axs[1].plot(s, y, "o")
    axs[1].grid()
    axs[1].legend(["аналитическое решение", "численное решение"])
    axs[1].set_ylabel("$y(s, t_1)$")
    axs[1].set_xlabel("$s$")

    print(np.max(abs(np.array(x_an_)-np.array(x))))
    print(np.max(abs(np.array(y_an_)-np.array(y))))
    plt.show()
