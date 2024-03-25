from abc import ABC, abstractmethod
from .problem import HypProblem
from .solve_method import SolveMethod


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
            #    f"\t ан.реш x:{x_an(self.s, self.t): .2f}, y:{y_an(self.s, self.t):.2f}"
               # f"({self.left.s:.2f}, {self.left.t:.2f}), ({self.right.s:.2f}, {self.right.t:.2f})"\
               # f"{type(self)}"
               # f"Node left: ({self.left.s: .2f}, {self.left.t: .2f})]"
               # f" [ Node right: ({self.right.s: .2f},{self.right.s: .2f})" \

    @abstractmethod
    def solve(self, problem: HypProblem, solve_method: SolveMethod):
        pass

    def get_inf(self):
        return [self.s, self.t, self.x, self.y]

    def get_old_point(self, problem: HypProblem, solve_method: SolveMethod):

        if not self.left.is_resolved:
            self.left.solve(problem, solve_method)

        sl, tl, xl, yl = self.left.get_inf()
        hl = self.t - tl

        if not self.right.is_resolved:
            self.right.solve(problem, solve_method)

        sr, tr, xr, yr = self.right.get_inf()
        hr = self.t - tr
        return (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr)


class NodeStart(Node):
    def create_parents_node(self, problem: HypProblem, solve_method: SolveMethod):
        if self.left is None:
            self.left = NodeStart(problem.C2*(problem.T0-self.t)+self.s, problem.T0)
            self.left.solve()

        if self.right is None:
            s =  problem.C1*(self.t - problem.T0)+self.s

            self.right = NodeStart(s,  problem.T0)
            if s > problem.S1:
                s1t0_node = NodeStart(problem.S1,  problem.T0)
                s1t0_node.solve(problem, solve_method)

                temp_node = NodeRight(self.s, self.t, self.left, s1t0_node)
                temp_node.solve(problem, solve_method)
                self.x = temp_node.x
                self.y = temp_node.y
                self.is_resolved = True
            else:
                self.right.solve(problem, solve_method)

    def solve(self, problem: HypProblem, solve_method: SolveMethod):
        if self.t == problem.T0:
            self.x = problem.x0(self.s)
            self.y = problem.y0(self.s)
        elif self.right is None or self.left is None:
            self.create_parents_node(problem, solve_method)
            if not self.is_resolved:
                (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = self.get_old_point(problem, solve_method)
                self.x, self.y = solve_method.start_solve(problem, self.s,  self.t, sl, tl, xl, yl, hl, sr, tr, xr, yr, hr)
        else:
            (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = self.get_old_point(problem, solve_method)
            self.x, self.y = solve_method.start_solve(problem, self.s, self.t, sl, tl, xl, yl, hl, sr, tr, xr, yr, hr)
        self.is_resolved = True



class NodeLeft(Node):
    def solve(self, problem: HypProblem, solve_method: SolveMethod):
        if self.t == problem.T0:
            self.x = problem.x0(self.s)
            self.y = problem.y0(self.s)
        else:
            (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = self.get_old_point(problem, solve_method)
            self.x, self.y = solve_method.left_solve(problem, self.s,  self.t, sl, tl, xl, yl, hl, sr, tr, xr, yr, hr)
        self.is_resolved = True


class NodeRight(Node):
    def solve(self, problem: HypProblem, solve_method: SolveMethod):
        if self.x is not None and self.y is not None:
            pass
        else:
            (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = self.get_old_point(problem, solve_method)
            self.x, self.y = solve_method.right_solve(problem, self.s,  self.t, sl, tl, xl, yl, hl, sr, tr, xr, yr, hr)

        self.is_resolved = True


class NodeFinish(Node):
    def solve(self, problem: HypProblem, solve_method: SolveMethod):
        (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = self.get_old_point(problem, solve_method)
        self.x, self.y = solve_method.finish_solve(problem, self.s,  self.t, sl, tl, xl, yl, hl, sr, tr, xr, yr, hr)
        self.is_resolved = True


class NoneCenter(Node):
    def solve(self, problem: HypProblem, solve_method: SolveMethod):
        (sl, tl, xl, yl, hl), (sr, tr, xr, yr, hr) = self.get_old_point(problem, solve_method)
        self.x, self.y = solve_method.center_solve(problem, self.s,  self.t, sl, tl, xl, yl, hl, sr, tr, xr, yr, hr)
        self.is_resolved = True