import matplotlib.pyplot as plt
import numpy as np
from .mesh import Mesh
from .nodes import Node

class Ploter:
    def plot_mesh(self, mesh: Mesh):
        start_nodes, left_nodes, right_nodes, center_nodes, finish_nodes = mesh.result

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

        plt.xlim([mesh.s0, mesh.s1])
        plt.ylim([mesh.t0, mesh.t1+0.5])
        plt.show()

    def plot_border_left(self, mesh: Mesh, x_an, y_an):
        start_nodes, left_nodes, right_nodes, center_nodes, finish_nodes = mesh.result
        fig, axs = plt.subplots(nrows=2, ncols=1)
        x = [node.x for node in left_nodes]
        y = [node.y for node in left_nodes]
        t = [node.t for node in left_nodes]
        axs[0].plot(t, x, "*")

        axs[0].plot(np.array(t), x_an(mesh.s0, np.array(t)))
        axs[1].plot(t, y, "*")
        axs[1].plot(np.array(t), y_an(mesh.s0, np.array(t)))

        plt.show()

    def plot_border_right(self, mesh: Mesh, x_an, y_an):
        start_nodes, left_nodes, right_nodes, center_nodes, finish_nodes = mesh.result
        fig, axs = plt.subplots(nrows=2, ncols=1)
        x = [node.x for node in right_nodes]
        y = [node.y for node in right_nodes]
        t = [node.t for node in right_nodes]
        axs[0].plot(t, x, "*")

        axs[0].plot(np.array(t), x_an(mesh.s1, np.array(t)))
        axs[1].plot(t, y, "*")
        axs[1].plot(np.array(t), y_an(mesh.s1, np.array(t)))

        plt.show()

    def plot_final(self, mesh: Mesh, x_an, y_an):
        s, x, y = mesh.get_XY_t1()
        x_an_ = [x_an(si, mesh.t1) for si in s]
        y_an_ = [y_an(si, mesh.t1) for si in s]

        fig, axs = plt.subplots(nrows=2, ncols=1)
        axs[0].set_title(f"Численное решение (узлов по s: {mesh.m})")
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

        axs[0].set_ylim([-0.1, 2])
        axs[1].set_ylim([-0.1, 2])

        print(np.max(abs(np.array(x_an_)-np.array(x))))
        print(np.max(abs(np.array(y_an_)-np.array(y))))
        plt.show()