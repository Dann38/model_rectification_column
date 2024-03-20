from node_sz import *

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
        self.right = MeshRight([self.h_down, self.h_up], self.s1, self.t1, inds[-1])
        self.center = MeshCenter(self.ds, [self.h_down, self.h_up], self.m, inds)
        self.finish = MeshFinish([self.h_down, self.h_up], inds)

        self.result = None

    def create_mesh_s(self):
        ds = (self.s1 - self.s0) / self.m
        h_down = ds / self.c2
        h_up = ds / self.c1
        return ds, h_down, h_up

    def create_mesh(self, revers_time=False):
        start_nodes, left_nodes, right_nodes, center_nodes, finish_nodes = [[] for i in range(5)]
        start_nodes = self.start.get_nodes()
        self.start.connect(start_nodes) # Хорошо работает

        left_nodes = self.left.get_nodes(start_nodes[0])
        self.left.connect(left_nodes)
        self.left.connect_start(left_nodes, start_nodes)

        right_nodes = self.right.get_nodes(start_nodes[-1], len(left_nodes))
        self.right.connect(right_nodes)
        self.right.connect_start(right_nodes, start_nodes)

        center_nodes = self.center.get_nodes(start_nodes[1:-1], len(left_nodes))
        self.center.connect(center_nodes)
        self.center.connect_start(center_nodes, start_nodes)
        self.center.connect_left(center_nodes, left_nodes)
        self.center.connect_right(center_nodes, right_nodes)

        self.left.connect_center(left_nodes, center_nodes)
        self.right.connect_center(right_nodes, center_nodes)

        finish_nodes = self.finish.get_nodes(left_nodes[-1], center_nodes[-1], right_nodes[-1])
        self.finish.connect(finish_nodes)

        self.finish.connect_center(finish_nodes, center_nodes)
        self.finish.connect_left(finish_nodes, left_nodes)
        self.finish.connect_right(finish_nodes, right_nodes)

        if revers_time:
            self.start.revers_time(start_nodes)
            self.left.revers_time()
        self.result = (start_nodes, left_nodes, right_nodes, center_nodes, finish_nodes)

    def solve(self):
        (start_nodes, left_nodes, right_nodes, center_nodes, finish_nodes) = self.result
        # print("получение стартовых узлов")
        for node in start_nodes:  # начальные условия
            node.solver()
        # print("получение оставшихся узлов")

        # print("решение узлов по уровням")
        for i, level_nodes in enumerate(center_nodes):
            left_nodes[i].solver()
            for node in level_nodes:
                node.solver()
            right_nodes[i].solver()

        for level in finish_nodes:
            for node in level:
                node.solver()

    def plot_mesh(self):
        start_nodes, left_nodes, right_nodes, center_nodes, finish_nodes = self.result

        def plot_arrow_node(node: Node):
            alpha = 0.8
            width_arrow = 0.01
            if node.left is not None:
                s = node.s
                t = node.t
                ds = (node.left.s - node.s)*alpha
                dt = (node.left.t - node.t)*alpha
                plt.arrow(s, t, ds, dt, width=width_arrow)
            if node.right is not None:
                s = node.s
                t = node.t
                ds = (node.right.s - node.s)*alpha
                dt = (node.right.t - node.t)*alpha
                plt.arrow(s, t, ds, dt, width=width_arrow)

        def plot_node(node: Node):
            color = ("r" if node.right is None else "b") if node.left is None else ("orange" if node.right is None else "g")
            plot_arrow_node(node)
            plt.scatter(node.s, node.t, color=color)

        for node in start_nodes:
            plot_node(node)
        for node in left_nodes:
            plot_node(node)
        for node in right_nodes:
            plot_node(node)
        for level in center_nodes:
            for node in level:
                plot_node(node)
        for level in finish_nodes:
            for node in level:
                plot_node(node)

        plt.xlim([self.s0, self.s1])
        plt.ylim([self.t0, self.t1+0.5])
        plt.show()

    def plot_border_left(self):
        start_nodes, left_nodes, right_nodes, center_nodes, finish_nodes = self.result
        fig, axs = plt.subplots(nrows=2, ncols=1)
        x = [node.x for node in left_nodes]
        y = [node.y for node in left_nodes]
        t = [node.t for node in left_nodes]
        axs[0].plot(t, x, "*")

        axs[0].plot(np.array(t), x_an(S0, np.array(t)))
        axs[1].plot(t, y, "*")
        axs[1].plot(np.array(t), y_an(S0, np.array(t)))

        plt.show()

    def plot_border_right(self):
        start_nodes, left_nodes, right_nodes, center_nodes, finish_nodes = self.result
        fig, axs = plt.subplots(nrows=2, ncols=1)
        x = [node.x for node in right_nodes]
        y = [node.y for node in right_nodes]
        t = [node.t for node in right_nodes]
        axs[0].plot(t, x, "*")

        axs[0].plot(np.array(t), x_an(S1, np.array(t)))
        axs[1].plot(t, y, "*")
        axs[1].plot(np.array(t), y_an(S1, np.array(t)))

        plt.show()

    def plot_final(self):
        s, x, y = self.get_XY_t1()
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

    def get_XY_t1(self):
        if self.result is None:
            self.solve()
        a1, a2, a3, a4, finish = self.result

        def approx(t_u, t_d, f_u, f_d):

            if t_d == t_u:
                return f_u
            return (f_d-f_u)*(T1-t_u)/(t_d-t_u) + f_u

        x_t1 = [approx(u.t, d.t, u.x, d.x) for u, d in zip(finish[0], finish[1])]
        y_t1 = [approx(u.t, d.t, u.y, d.y) for u, d in zip(finish[0], finish[1])]
        s_ = [n.s for n in finish[0]]

        return s_, x_t1, y_t1


class MeshStart:
    def __init__(self, ds, h,  m, s0, t0):
        self.ds = ds
        self.h = h
        self.m = m
        self.s0 = s0
        self.t0 = t0
        self.inds = self.neighbor_is_left()

    def get_nodes(self):
        nodes = []
        t_ = self.t0
        s_ = self.s0

        for i in range(self.m+1):
            nodes.append(NodeStart(s_, t_))
            t_ = (t_ + self.h[0]) % (self.h[0]+self.h[1])
            s_ += self.ds

        return nodes

    def connect(self, nodes):

        for i in range(1, len(nodes)):
            nodes[i].left = nodes[i-1] if self.inds[i-1] == 1. else None
            nodes[i-1].right = nodes[i] if self.inds[i-1] == 0. else None

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

    def get_nodes(self, start_nodes, count_level):
        nodes = [[NoneCenter(node.s, node.t+self.dt*i) for node in start_nodes] for i in range(1, count_level+1)]
        return nodes

    def connect(self, nodes):
        count_level = len(nodes)
        for level in range(1, count_level):
            for i, ind in enumerate(self.inds[1:-1]):
                nodes[level][i + 1].left = nodes[level][i] if ind == 1. else nodes[level - 1][i]
                nodes[level][i].right = nodes[level - 1][i + 1] if ind == 1. else nodes[level][i + 1]

        for i, ind in enumerate(self.inds[1:-1]):
            if ind == 1.:
                nodes[0][i + 1].left = nodes[0][i]
            else:
                nodes[0][i].right = nodes[0][i + 1]

    def connect_start(self, nodes, start_nodes):
        for i, ind in enumerate(self.inds[1:-1]):
            if ind == 1.:
                nodes[0][i].right = start_nodes[i + 2]
            else:
                nodes[0][i + 1].left = start_nodes[i+1]
        if self.inds[-1] == 1.:
            nodes[0][-1].right = start_nodes[-1]

    def connect_left(self, nodes, left_nodes):
        # Стыковка левых точек ----------------------------------------------------------------------------------------
        for i in range(len(nodes)):
            nodes[i][0].left = left_nodes[i]

    def connect_right(self, nodes, right_nodes):
        # Стыковка правых точек ---------------------------------------------------------------------------------------
        if self.inds[-1] == 1.:
            for i in range(1, len(nodes)):
                nodes[i][-1].right = right_nodes[i-1]
        else:
            for i in range(len(nodes)):
                nodes[i][-1].right = right_nodes[i]


class MeshLeft:
    def __init__(self, h, s0, t1):
        self.h = h
        self.s0 = s0
        self.t1 = t1

    def get_nodes(self, start_node):
        dt = self.h[0]+self.h[1]
        t = np.arange(start_node.t+dt, self.t1-dt, dt)
        nodes = [NodeLeft(self.s0, ti) for ti in t]
        return nodes

    def connect(self, nodes):
        for i, node in enumerate(nodes[1:]):
            node.left = nodes[i]

    def connect_start(self, nodes, start_nodes):
        nodes[0].left = start_nodes[0]
        nodes[0].right = start_nodes[1]

    def connect_center(self, nodes, center_nodes):
        for i in range(len(nodes)-1):
            nodes[i+1].right = center_nodes[i][0]


class MeshRight:
    def __init__(self, h, s1, t1, ind):
        self.h = h
        self.s1 = s1
        self.t1 = t1
        self.ind = ind

    def get_nodes(self, start_node, count_nodes):
        dt = self.h[0] + self.h[1]
        t = [ start_node.t+i*dt for i in range(1, count_nodes+1)]
        nodes = [NodeRight(self.s1, ti) for ti in t]
        return nodes

    def connect(self, nodes):
        for i, node in enumerate(nodes[1:]):
            node.right = nodes[i]

    def connect_start(self, nodes, start_nodes):
        nodes[0].right = start_nodes[-1]
        if self.ind == 0:
            nodes[0].left = start_nodes[-2]

    def connect_center(self, nodes, center_nodes):
        if self.ind == 1.:
            for i in range(len(nodes)):
                nodes[i].left = center_nodes[i][-1]
        else:
            for i in range(1, len(nodes)):
                nodes[i].left = center_nodes[i-1][-1]


class MeshFinish:
    def __init__(self, h, inds):
        self.dt = h[1]+h[0]
        self.inds = inds

    def get_nodes(self, left_node, center_nodes, right_node):
        finish_nodes = [[],[]]
        for i in [0, 1]:
            dt = (1+i)*self.dt
            finish_nodes[i] = ([NodeLeft(left_node.s, left_node.t+dt)] +
                        [NoneCenter(node.s, node.t+dt) for node in center_nodes] +
                        [NodeRight(right_node.s, right_node.t+dt)])


        return finish_nodes

    def connect(self, nodes):
        nodes[1][1].left = nodes[1][0]
        nodes[1][0].right = nodes[0][1]
        nodes[1][0].left = nodes[0][0]
        nodes[0][1].left = nodes[0][0]
        for i, ind in enumerate(self.inds[1:-1]):
            if ind == 1.:
                nodes[0][i + 2].left = nodes[0][i + 1]
                nodes[1][i + 2].left = nodes[1][i + 1]
                nodes[1][i + 1].right = nodes[0][i + 2]  # соединение с нижним
            else:
                nodes[0][i + 1].right = nodes[0][i + 2]
                nodes[1][i + 1].right = nodes[1][i + 2]
                nodes[1][i + 2].left = nodes[0][i + 1]  # соединение с нижним
        if self.inds[-1] == 1.:
            nodes[0][-1].left = nodes[0][-2]
            nodes[1][-1].left = nodes[1][-2]
            nodes[1][-2].right = nodes[0][-1]
        else:
            nodes[0][-2].right = nodes[0][-1]
            nodes[1][-2].right = nodes[1][-1]
            nodes[1][-1].left = nodes[0][-2]
        nodes[1][-1].right = nodes[0][-1]

    def connect_center(self, nodes, center_nodes):
        nodes[0][0].right = center_nodes[-1][0]
        for i, ind in enumerate(self.inds[1:-1]):
            if ind == 1.:
                nodes[0][i + 1].right = center_nodes[-1][i + 1]  # соединение с нижним
            else:
                nodes[0][i + 2].left = center_nodes[-1][i]  # соединение с нижним
        if self.inds[-1] == 0.:
            nodes[0][-1].left = center_nodes[-1][-1]

    def connect_left(self, nodes, left_nodes):
        nodes[0][0].left = left_nodes[-1]

    def connect_right(self, nodes, right_nodes):
        if self.inds[-1] == 1.:
            nodes[0][-2].right = right_nodes[-1]
        nodes[0][-1].right = right_nodes[-1]

